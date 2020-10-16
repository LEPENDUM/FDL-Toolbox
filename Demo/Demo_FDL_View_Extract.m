%-------------------------------------------------------------------------%
%   Demo script for the FDL completion/super-resolution/demosaicing
%   algorithm to extract a Light Field from Lytro camera RAW data.
%   See the paper: M. Le Pendu and A. Smolic, "High Resolution Light Field
%   Recovery with Fourier Disparity Layer Completion, Demosaicing, and
%   Super-Resolution", ICCP 2020.
%
%   The Lytro RAW data must be converted first to generate incomplete views
%   using an external toolbox (see instructions below).
%   The scripts generates the final central view and focal stack.
%
%   Instructions:
%   -1. Capture Lytro image or download RAW Lytro images from existing
%   datasets (e.g.: <a href="matlab: web('http://clim.inria.fr/research/LowRank2/datasets/datasets.html')">INRIA Dataset</a>, <a href="matlab: web('https://jpeg.org/plenodb/lf/epfl/')">EPFL Dataset</a>).
%
%   -2. Convert RAW data (.LFR file) into incomplete views with the Light
%   Field toolbox available here: 
%           https://github.com/V-Sense/LFToolbox-CLIM_VSENSE
%   The decoding options should be set as:
%   DecodeOptions.OptionalTasks = 'ColourCorrect';
%   DecodeOptions.ResampMethod = 'none';
%   The other options should be left as default.
%   Settting the ResampMethod to 'none' generates images in the correct
%   format (e.g. avoid any destructive operations such as resampling
%   bilinear/bicubic interpolations, generate masks indicating the missing
%   pixels, ...).
%
%   -3. Place the decoded mat file (finishing with '__Decoded.mat') in the
%   Demo directory.
%
%   -4. Set the variable LFName and run the script.
%
%-------------------------------------------------------------------------%

LFDir = fullfile(fileparts(mfilename('fullpath')), '/');
useGPU = parallel.gpu.GPUDevice.isAvailable;
saveFDL = true;	%Save the resulting FDL
reload = false;	%Load an existing FDL model and skip processing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Input Light Field parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LFName = 'Bee_2';           DispRange = [];
%LFName = 'Danger_de_Mort';  DispRange = [];
%LFName = 'Vespa';           DispRange = [];
%LFName = 'Fruits';          DispRange = [1 5];
%LFName = 'Bumblebee';       DispRange = [-10 10];
%LFName = 'Field';           DispRange = [];

LFnameSuffix='__Decoded';

if(~reload)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parameters for FDL calibration and construction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambdaCalib=1;          % Regularization parameter for FDL calibration.
numIter = 30;           % Number of iteration of the FDL completion algorithm.
lambdaCstrBase=0.03;    % l2 regularization parameter.
lambdaColor = 1;        % Color regularization parameter.
sigSpatial=.4;          % Spatial blur std. dev. (w.r.t low resolution image scale).
sigAngular=1.3;         % Angular blu std. dev. (assuming normalized baseline).
scaleX = 2;             % Super-resolution horizontal scale
scaleY = 2;             % Super-resolution vertical scale
useAngularPrecond=true;

sigX=sigSpatial;
sigY=sigSpatial*2/sqrt(3);%Extracted views are stretched vertically by a factor 2/sqrt(3) due to the hexagonal lenslet sampling.

padSizeX=15;			% Number of padded pixels on left and right borders.
padSizeY=15;			% Number of padded pixels on top and bottom borders.
gammaOffset = 0;		% Offset applied before inverse gamma correction.
MaskSkip=1;             % Number of pixels to discard on teach border of the mask.



%% Load Light Field data
tic;fprintf(['Loading input views for Light Field ''' LFName ''' ...']);
load([LFDir LFName LFnameSuffix '.mat']);
LF=single(permute(LF,[3 4 5 1 2]))/65535;
ViewWeights = DecodeOptions.WIAvgPerView / max(DecodeOptions.WIAvgPerView(:));
ColMatrix = DecodeOptions.ColourMatrix;
numInitViews = size(LF,4);
imgSize = [size(LF,1),size(LF,2)];
nChan = DecodeOptions.NColChans;
nWChan = DecodeOptions.NWeightChans;

t=toc;fprintf(['\b\b\b (done in ' num2str(t) 's)\n']);


%% Prepare sampling variables
oddY = mod(imgSize(1),2);%used to foce an even number of horizontal lines (for accurate handling of hexagonal sampling with FDL 2x super-resolution).
fullResX = padSizeX*2+imgSize(2);
fullResY = padSizeY*2+imgSize(1) + oddY;
hexaSampling = DecodeOptions.Sampling.HexShiftStart;%1->odd rows shifted by 1/2 pixel / 2->even rows shifted by 1/2 pixel
%frequency grid:
[wx,wy,xC,yC] = GenerateFrequencyGrid(fullResX,fullResY);


%% Pre-process input data (Padding/Windowing + Fourier Tramsform )
tic;fprintf('Pre-process RAW data: Padding/Windowing/Fourier Transform)...');
% Prepare mask data with zero padding:
bp = BorderParams('X',padSizeX+MaskSkip,'T',padSizeY+MaskSkip,'B',padSizeY+oddY+MaskSkip,'padVal',0,'window','none');
Mask = padImgs(logical(LF(1+MaskSkip:end-MaskSkip,1+MaskSkip:end-MaskSkip,nChan+1:nChan+nWChan,:)),bp,'disableGPU',~useGPU);

% Compute Fourier transfom of the input views for calibration (including padding and windowing):
LF = single(LF(:,:,1:nChan,:));
bp = BorderParams('X',padSizeX,'T',padSizeY,'B',padSizeY+oddY,'padVal','replicate','window','hann');
ViewsFFT = fftImgs(LF,bp,'dimsOrderOut',[4 3 1 2],'hexaSampling',hexaSampling,'disableGPU',~useGPU);

% Prepare input views with padding:
bp.window='none';
LF = padImgs(single(LF(:,:,1:nChan,:)),bp,'disableGPU',~useGPU);

%Swap the hexagonal sampling type in the case of an odd padding so that hexaSampling
%defines the parity of row indices for the images with padded borders.
hexaSampling = 2 - mod(hexaSampling + mod(bp.T,2), 2);%swap values 1 and 2 if bp.T is odd.


t=toc;fprintf(['\b\b\b (done in ' num2str(t) 's)\n']);

%% FDL Calibration
numLayersInit=20;
numLayers=numLayersInit;

fprintf('FDL Calibration...');
tic;[~,~,U0,V0,D0]=CalibrateFDL_UVD_gpu(ViewsFFT, wx, wy, numLayersInit,lambdaCalib);toc;
%tic;[~,~,Urgb,Vrgb,Drgb,~,~,~,~,~,residN2,gradN2]=CalibrateFDL_UVD_gpu(ViewsFFT, wx, wy, numLayers,lambdaCalib,U0,V0,D0,'SeparateChannels',true);toc;
tic;[~,~,Urgb,Vrgb,Drgb,~,~,~,~,~,residN2,gradN2]=CalibrateFDL_UVD_gpu(ViewsFFT, wx, wy, numLayers,lambdaCalib,U0,V0,D0,'SeparateChannels',true,'UseSignFlip',false,'numIterLowFreq',0);toc;
%% Select Disparity values within the specified disparity Range.
if(~isempty(DispRange))
    minD = min(Drgb,[],2);
    maxD = max(Drgb,[],2);
    if(DispRange(1)<max(minD(:))), DispRange(1) = max(minD(:));end
    if(DispRange(2)>min(maxD(:))), DispRange(2) = min(maxD(:));end
    
    Id_avg = interp1(mean(Drgb,3), 1:numLayersInit, linspace(DispRange(1), DispRange(2), numLayersInit), 'pchip');
    minId = min(Id_avg(:));
    maxId = max(Id_avg(:));
    numLayers = 1+ceil(maxId-minId);
    clear Id_avg
    Drgb_ = zeros(1,numLayers,nChan);
    for ch=1:nChan
        Drgb_(:,:,ch) = interp1(Drgb(:,:,ch), linspace(minId, maxId, numLayers), 'pchip');
    end
    Drgb=Drgb_;
    clear Drgb_
end

lambdaConstruct = lambdaCstrBase * numLayersInit / numLayers ;


%% Scale angular coordinates according to sensor sampling.
UnitRefU=DecodeOptions.Sampling.ST(:,1)/DecodeOptions.Sampling.SubPixAccuracy;
UnitRefV=DecodeOptions.Sampling.ST(:,2)/DecodeOptions.Sampling.SubPixAccuracy;

[deltaU,deltaV]=approxRegGrid(Urgb(:,:,2),Vrgb(:,:,2),UnitRefU,UnitRefV);
deltaU = abs(deltaU);
deltaV = abs(deltaV);

%For test without reduction of chromatic aberration.
%Urgb=mean(Urgb,3); Vrgb=mean(Vrgb,3); Drgb=mean(Drgb,3);

%[FDL,vidF]=FDL_Complete_SR(permute(LF,[4 3 1 2]), permute(Mask,[4 3 1 2]), U, V, D, useAngularPrecond, deltaU, deltaV, [scaleX scaleY], lambdaConstruct*5/6, lambdaColor, [], [], [sigX, sigY, sigAngular], numIter, hexaSampling, ViewWeights, {}, 'HR-FDL');%for precond1 (upsampling first)
%%{
[FDL,vidF]=FDL_Complete_SuperRes(LF, Mask, Urgb, Vrgb, Drgb, lambdaConstruct, [scaleX scaleY], [sigX, sigY],...
    'useAngPrecond',useAngularPrecond,...
    'deltaUV',[deltaU deltaV],...
    'lambdaColor',lambdaColor,...
    'sigAngular',sigAngular,...
    'numIter',numIter,...
    'hexaSampling',hexaSampling,...
    'ViewWeights',ViewWeights,...
    'outputMode','HR-FDL');
%}

%% Resize final results for correct aspect/ratio
XExtRatio=2/sqrt(3);
sizeFDLX = ceil(size(FDL,2)*XExtRatio);
sizeFDLY = size(FDL,1);
FDL=XExtRatio*padarray(FDL,[0 floor((sizeFDLX-size(FDL,2))/2) 0 0],0);
FDL(:,end+mod(sizeFDLX,2),:,:)=0;

padSizeXExt = [floor(padSizeX*scaleX*XExtRatio) ceil(padSizeX*scaleX*XExtRatio)];
padSizeYExt = scaleY*padSizeY;
U = scaleX*XExtRatio*mean(Urgb,3);
V = scaleY*mean(Vrgb,3);
Disps = mean(Drgb,3);

%% Create Render Model object
%RenderAppMain(FDL, [sizeFDLY sizeFDLX], [padSizeXExt(1) padSizeXExt(2) padSizeYExt padSizeYExt], mean(D,3), scaleX*XExtRatio*mean(U,3), scaleY*mean(V,3),0);
RMod = RenderModel(FDL, [sizeFDLY sizeFDLX], [padSizeXExt(1) padSizeXExt(2) padSizeYExt padSizeYExt], Disps);

if(saveFDL)
    RMod.saveFDL([LFDir LFName '.fdl'], U, V);
    save([LFDir LFName '.fdl'], 'ColMatrix', '-append');
end

doPostProc=true;

else%load existing FDL
    load([LFDir LFName '.fdl'], '-mat');
    doPostProc = ~isempty(who('-file', [LFDir LFName '.fdl'], 'ColMatrix'));
    nChan = size(FDL,3);
    RMod = RenderModel(FDL, fullSize, crop, Disps, 0, [], [], useGPU);
end

%% Render final images (central view and Focal Stack with full aperture) and post-process colors (color matrix transform in linearized RGB space).
%Render Central View:
RMod.setPosition(0,0);
RMod.renderImage();
imgDimsFinal = size(RMod.Image);
CentralView = RMod.getInternalImage(true);
if(doPostProc)
    CentralView = reshape(sRGB_gammaEncode(reshape(sRGB_gammaDecode(CentralView),[],nChan)*ColMatrix), imgDimsFinal);
end
figure,imshow(CentralView)

%Render Refocused images
numRefocus = numel(Disps);
FocalStack = zeros([imgDimsFinal numRefocus]);
RMod.setApShape('disk');
RMod.computeAperture;
RMod.setRadius(sqrt(max(U(:).^2+V(:).^2)));
for idRefocus=1:numRefocus
    RMod.setFocus(Disps(idRefocus));
    RMod.renderImage();
    FocalStack(:,:,:,idRefocus) = RMod.getInternalImage(true);
    if(doPostProc)
        FocalStack(:,:,:,idRefocus) = reshape(sRGB_gammaEncode(reshape(sRGB_gammaDecode(FocalStack(:,:,:,idRefocus)),[],nChan)*ColMatrix), imgDimsFinal);
    end
end
implay(FocalStack);