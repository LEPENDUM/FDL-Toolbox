%-------------------------------------------------------------------------%
%          Demo script for the FDL super-resolution algorithm
%
%   See the paper: M. Le Pendu and A. Smolic, "High Resolution Light Field
%   Recovery with Fourier Disparity Layer Completion, Demosaicing, and
%   Super-Resolution", ICCP 2020.
%-------------------------------------------------------------------------%

LFDir = [fileparts(mfilename('fullpath')) '/'];
useGPU = parallel.gpu.GPUDevice.isAvailable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Input Light Field parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LFName = 'Illum_Field';
uRange = [3:11];		%list of u indices of the images to load from the files.
vRange = [3:11];		%list of v indices of the images to load from the files.
crop = [1 1 1 1];		%croped pixels form the input images [left,right,top,bottom].
padType = 'symmetric';	%'replicate' or 'symmetric'
Windowing = true;		%Use window to decrease intensity of padded pixels.
padSizeX=15;			%Number of padded pixels on left and right borders.
padSizeY=15;			%Number of padded pixels on top and bottom borders.
linearizeInput = true;	%Use inverse gamma correction to process images in linear space.
gammaOffset = 0;		%Offset applied before inverse gamma correction.


nU = length(uRange);
nV = length(vRange);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parameters for FDL calibration + construction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambdaCalib = 1;
lambdaConstruct_L2 = 4;    %higher -> more reduction of noise/artifacts.
lambdaConstruct_Color = 1; %0: no color regularisation (faster) / >0: reduce color noise/artifacts
SuperRes_Factor = 2;
SuperRes_SpatSigma = 0.33;
numLayers = 20;

%% Load Light Field data
tic;fprintf('Loading input views...');
LF = loadLF([LFDir LFName], '', 'png', uRange, vRange, crop);
LF = double(reshape(LF,size(LF,1),size(LF,2),size(LF,3),[]))/255;
imgSize = [size(LF,1),size(LF,2)];
nChan = size(LF,3);
if(linearizeInput)
    LF = BT709_gammaDecode(LF+gammaOffset);
end
t=toc;fprintf(['\b\b\b (done in ' num2str(t) 's)\n']);


%% Prepare data
fullResX = padSizeX*2+imgSize(2);
fullResY = padSizeY*2+imgSize(1);
%frequency grid:
[wx,wy,xC,yC] = GenerateFrequencyGrid(fullResX,fullResY);
%Border parameters for padding/windowing of padded borders:
bp = BorderParams('X',padSizeX,'Y',padSizeY,'padVal',padType);
if(~Windowing)
    bp.window='none';
end

%Fourier Transform (with padding/windowing) of input views.
tic;fprintf('Computing Fourier Transform of input views...');
ViewsFFT = fftImgs(LF,bp,'dimsOrderOut',[4 3 1 2],'disableGPU',~useGPU);
t=toc;fprintf(['\b\b\b (done in ' num2str(t) 's)\n']);


%% FDL Calibration
tic;fprintf('FDL Calibration...');
if(useGPU)
    [Px,Py,U,V,D,~,~,~,~,~,residN2,gradN2]=CalibrateFDL_UVD_gpu(ViewsFFT, wx, wy, numLayers, lambdaCalib);
else
    [Px,Py,U,V,D,~,~,~,~,~,residN2,gradN2]=CalibrateFDL_UVD_cpu(ViewsFFT, wx, wy, numLayers, lambdaCalib);
end
%Find vertical and horizontal baseline (deltaU, deltaV) by fitting the U,V positions to an integer grid of views.
[uGridReg,vGridReg] = meshgrid(uRange-mean(uRange), vRange-mean(vRange));
[deltaU,deltaV,u0,v0]=approxRegGrid(U,V,uGridReg,vGridReg);
%Scaling and centering of the parameters to match approximately the integer grid of views with a baseline of 1 (this step has no effect on the final images).
divUV = (abs(deltaU)+abs(deltaV))/2;
U = (U - u0)/ divUV; u0=0; deltaU = deltaU / divUV;
V = (V - v0)/ divUV; v0=0; deltaV = deltaV / divUV;
D = D * divUV;

t=toc;fprintf(['\b(' num2str(t) 's)\n']);

%% FDL Construction
tic;fprintf('FDL Construction...');
if(useGPU)
    FDL = ComputeFDL_SuperRes_gpu(ViewsFFT, U, V, D, lambdaConstruct_L2, SuperRes_Factor, SuperRes_SpatSigma,...
        'lambdaColor', lambdaConstruct_Color, ...
        'deltaUV',[deltaU -deltaV]);
else
    error('FDL Super resolution not implmented for CPU.');
end
t=toc;fprintf(['(' num2str(t) 's)\n']);

%% Start Rendering Application
bp_HR=BorderParams('X',padSizeX*SuperRes_Factor,'Y',padSizeY*SuperRes_Factor);
fullRes_HR = SuperRes_Factor * [fullResY fullResX];
D_HR = D * SuperRes_Factor;
RenderAppMain(FDL, fullRes_HR, bp_HR, D_HR,U,V, [],{linearizeInput,gammaOffset}, useGPU);
