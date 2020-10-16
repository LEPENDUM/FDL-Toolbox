%-------------------------------------------------------------------------%
%   Demo script for the FDL construction algorithm (no super-resolution)
%
%   See the paper: M. Le Pendu, C. Guillemot and A. Smolic, "A Fourier
%   Disparity Layer Representation for Light Fields", Transactions on Image
%   Processing, 2019.
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
lambdaConstruct_L2 = 1;
lambdaConstruct_2ndOrderView = 0.1;
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
[U,V,D]=ScaleParams(U,V,D,nU,nV); %Scaling and centering of the parameters to match approximately an integer grid of views (this step has no effect on the final images).
t=toc;fprintf(['\b(' num2str(t) 's)\n']);

%% FDL Construction
tic;fprintf('FDL Construction...');
if(useGPU)
    FDL = ComputeFDL_gpu(ViewsFFT, wx, wy, U, V, D, [lambdaConstruct_L2 lambdaConstruct_2ndOrderView]);
else
    FDL = ComputeFDL_cpu(ViewsFFT, wx, wy, U, V, D, [lambdaConstruct_L2 lambdaConstruct_2ndOrderView]);
end
t=toc;fprintf(['(' num2str(t) 's)\n']);

%% Start Rendering Application
RenderAppMain(FDL, [fullResY fullResX], bp, D,U,V, [],{linearizeInput,gammaOffset}, useGPU);
