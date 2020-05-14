%Iterative FDL completion function.
%The completion can be performed jointly with super-resolution/demosaicing/deconvolution and with different FDL parameters per channel for chromatic aberations.
%
%---------------------------  Inputs arguments ----------------------------
% - Views:
%       4D Array of incomplete and low resolution Light Field views in the pixel domain
%       with dimensions: 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color channels, 4.Views.
%
% - Mask:
%       Binary Mask of known (true) and unknown (false) pixels.
%       The spatial resolution (first two dimensions) must be the same as in the 'Views' array.
%       The other dimensions (Views, Color channels) can be either the same as
%       the 'Views' array (e.g. different mask for each view and/or channel).
%       or singleton dimensions (same mask for all the views and/or channel).
%
%  -U, V:
%       Arrays of horizontal (U) and vertical (V) view positions (see FDL calibration tools to determine U,V automatically).
%       The dimensions of U,V can be (#Views x 1) or (#Views x #channels) to use different view positions for each color channel.
%
%  -D:
%       Array of disparity values. (see FDL calibration tools to determine U,V automatically).
%       The dimensions of D can be (1 x #layers) or (1 x #layers x #channels) to use different view positions for each color channel.
%
%  -lambda:
%       l2 regularization parameter.
%
%  -SRFact (default=1):
%       Super-resolution factor: SRFact(1) = horizontal factor.
%                                SRFact(2) = vertical factor (if not defined, same as horizontal factor).
%
%  -dwnFilt (default=0):
%      Spatial filter for deconvolution. It must be specified as follows:
%       -For a continuous gaussian filter with std. dev. parameter sigma (relative to the input image resolution) use either:
%           -A numeric value dwnFilt=sigma (or vector dwnFilt=[sigma_x,sigma_y] for different std. dev. horizontally and vertically).
%           -Or, a cell array with dwnFilt{1}='exactGauss', dwnFilt{2}=sigma (or dwnFilt{2}=[sigma_x,sigma_y]).
%       -For a discrete box filter of size (SRFactX x SRFactY), use a cell array with dwnFilt{1}='box'.
%       -For a custom discrete filter defined according to the output spatial resolution, use a cell array with:
%           -dwnFilt{1}='custom'.
%           -dwnFilt{2}=dwnFiltKernel -> The 2D kernel defining the filter.
%           -dwnFilt{3}=dwnFiltKernelOX -> Horizontal index of the kernel's origin (with indices from 1 to size(dwnFiltKernel,2)).
%           -dwnFilt{4}=dwnFiltKernelOY -> Vertical index of the kernel's origin (with indices from 1 to size(dwnFiltKernel,1)).
%           If dwnFilt{3} and dwnFilt{4} are not given, the origin is set to the center by default.
%
%  -varargin: Additional optional arguments given as (Name,Value) pairs:
%       -numIterMax (default=30):
%           Number of iterations.
%       -Px: Matrix of parameters for the horizontal dimension (default = [], i.e. automatically set as the outer product of U and D).
%       -Py: Matrix of parameters for the vertical dimension (default = [], i.e. automatically set as the outer product of V and D).
%            Px and Py may also include a 3rd dimension for color components.
%       -lambdaColor (default=1):
%           Color Regularization parameter.
%       -deltaUV (default=1):
%           Light Field baseline, i.e. difference in U,V coordinates between consecutive horizontal views. Required for angular filter and angular pre-conditionning.
%           'deltaUV' can also be a vector [deltaU, deltaV] to specify different vertical and horizontal baselines.
%       -sigAngular (default=0):
%           Standard deviation parameter for gaussian anglar blur model.
%       -useAngPrecond (default=true):
%           Set to true to use angular preconditioning. The amount of angular pre-conditioning is contoled with deltaUV parameters.
%       -hexaSampling: Value indicating the sampling format of the input images:
%           0 (default) -> Square input sampling.
%           1 -> Hexagonal input sampling with odd horizontal lines shifted to the right (unknown top left corner pixel).
%           2 -> Hexagonal input sampling with even horizontal lines shifted to the right (known top left corner pixel).
%           For Hexagonal samplings, the vertical resolution must be a multiple of 2 and the horizontal scale must be either 1 or a multiple of 2.
%       -ViewWeights:
%           Vector of weights for each view (default=[], i.e. use same weight=1 for each view).
%       -SRFactPrecond:
%           Preconditioning scale factor (default=[], i.e. automatically set equal to the SRFact).
%           This parameter may be set different from SRFact
%           e.g. set to 1 to disable spatial preconditioning,
%                set higher than SRFact to use preconditioning on already upsampled input.
%           'SRFactPrecond' can also be a vector [SRFactPrecondX SRFactPrecondY] to specify different values horizontally and vertically.
%       -nonNegative (default=true):
%           If set to true, prevents negative values at each iteration of the FDL construction.
%       -freqBlockSz (default=1024):
%           Number of frequencies processed in parallel.
%           It may be adjusted to prevent 'out of memory' depending on available GPU memory, number of views, and number of layers.
%       -nbImgGPU: Number of Images treated in parallel for GPU computation (default=50, may need adjustements depending on gpu memory available).
%       -seqIterStep (default=0):
%           -If seqIterStep > 0:  Generate sequence of images showing result of one view (with central index) at every 'seqIterStep' iteration.
%           -If seqIterStep <= 0: Skip image sequence generation (default option).
%       -ViewsRef (default=[]):
%           Reference Light field used to compute mean square error at each iteration (MSE is not computed if ViewsRef is empty or has wrong dimensions).
%       -outputMode:
%           -'LR': return the low resolution light field resulting from super-resolution and downsampling (i.e. low resolution FDL filtering preserving spatial aliasing).
%           -'HR': return the reconstructed High resolution light field (no residual back-projection).
%           -'HR-FDL' (default): return the high resolution FDL model.
%
%-------------------------------- Outputs ---------------------------------
% - F: Output in the Fourier domain as specified by the 'outputMode' argument.
%      For a FDL output ('HR-FDL' mode), the dimensions are in the order : 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color channels, 4.Layers.
%      For views output, the dimensions are in the order: 1.Views, 2.Color channels, 3.Vertical axis (X), 4.Horizontal axis (X).
%
% - IterInfo:
%     Structure containing per-iteration data:
%       -IterInfo.seqIter: Video showing result of one view (with central index) at every 'seqIterStep' iteration.
%       -IterInfo.seqIterStep: Step of the video in number of iterations.
%       -IterInfo.MSE: Mean square error w.r.t. the Reference Light field ViewsRef (if ViewsRef is specified as optional input argument and has the same size as Views).
%

function [F,IterInfo] = FDL_Complete_SuperRes(Views, Mask, U, V, D, lambda, SRFact, dwnFilt, varargin)

p = inputParser;

addParameter (p,'outputMode',       'HR-FDL',   @(x)any(validatestring(x,{'LR','HR','HR-FDL'})));
addParameter (p,'lambdaColor',      1,          @isnumeric);
addParameter (p,'hexaSampling',     0,          @(x)isnumeric(x)&&(x==0||x==1||x==2));
addParameter (p,'deltaUV',          1,          @(x)isnumeric(x) && numel(x)>=1 && numel(x)<=2);
addParameter (p,'useAngPrecond',    true,       @islogical);
addParameter (p,'sigAngular',       0,          @isnumeric);
addParameter (p,'ViewWeights',      [],         @isnumeric);
addParameter (p,'SRFactPrecond',    [],         @(x)isnumeric(x) && numel(x)<=2);
addParameter (p,'freqBlockSz',      1024,       @isnumeric);
addParameter (p,'Px',               [],         @isnumeric);
addParameter (p,'Py',               [],         @isnumeric);

addParameter (p,'numIterMax',       30,         @(x)isnumeric(x)&&x>0);
addParameter (p,'nonNegative',      true,       @islogical);
addParameter (p,'disableGPU',       false,      @islogical);
addParameter (p,'nbImgGPU',         50,         @(x)isnumeric(x)&&x>0);

addParameter (p,'seqIterStep',      0,          @isnumeric);
addParameter (p,'ViewsRef',         [],         @isnumeric);

parse(p,varargin{:});
outputMode = p.Results.outputMode;
lambdaColor = p.Results.lambdaColor;
hexaSampling = p.Results.hexaSampling;
deltaUV = p.Results.deltaUV;
useAngPrecond = p.Results.useAngPrecond;
sigAng = p.Results.sigAngular;
ViewWeights = p.Results.ViewWeights;
SRFactPrecond = p.Results.SRFactPrecond;
freqBlockSz = p.Results.freqBlockSz;
Px = p.Results.Px;
Py = p.Results.Py;

numIterMax = p.Results.numIterMax;
nonNegative = p.Results.nonNegative;
disableGPU = p.Results.disableGPU;
nbImgGPU = p.Results.nbImgGPU;
seqIterStep = p.Results.seqIterStep;
ViewsRef = p.Results.ViewsRef;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useGPU = parallel.gpu.GPUDevice.isAvailable && ~disableGPU;

if(~exist('SRFact','var') || isempty(SRFact))
    SRFact = 1;
end
if(~exist('dwnFilt','var') || isempty(dwnFilt))
    dwnFilt = 0;
end

%Size Variables
imgSize = [size(Views,1),size(Views,2)];
nChan = size(Views,3);
numImages = size(Views,4);
if(size(Mask,1)~=imgSize(1) || size(Mask,2)~=imgSize(2))
    error('The spatial resolution of the Mask (first two dimensions) should be equal to that of the input Views.');
end
if(size(Mask,3)~=nChan && size(Mask,3)~=1)
    error('The number of channels of the Mask (3rd dimension) should be either 1 or equal to that of the input Views.');
end
if(size(Mask,4)~=numImages && size(Mask,4)~=1)
    error('The number of views of the Mask (4th dimension) should be either 1 or equal to that of the input Views.');
end
%Initialization
F = single(Views);
Views = bsxfun(@times, Mask, Views); %The masked views will be used to reset the known pixels to their original values at each iteration.

%Prepare iterInfo data:
IterInfo=[];
computeMSE = isequal(size(ViewsRef),size(Views));
if(~isempty(ViewsRef) && ~computeMSE)
    warning('The reference light field ViewsRef specified for MSE computation should be the same size as the input light field -> MSE won''t be computed.')
end
if(computeMSE)
    IterInfo.MSE = zeros(1,numIterMax);
end

if(seqIterStep>0)
    seqIterStep = round(seqIterStep);
    IterInfo.seqIterStep = seqIterStep;
    IterInfo.seqIter = zeros(imgSize(1),imgSize(2),nChan,ceil((numIterMax-1)/seqIterStep));
    computeVideo = true;
else
    computeVideo = false;
end

for it=1:numIterMax

    if(it==numIterMax)
        outputModeCur = outputMode;
    else
        outputModeCur = 'LR';
    end
    
%% Forward Fourier transform and Hexagonal to Square resampling:
    fprintf(['iter#' num2str(it) ': ']);
    tic;fprintf('Forward FFT...');
    F = fftImgs(F, [], 'dimsOrderOut', [4 3 1 2], 'hexaSampling',hexaSampling,'nbImgGPU',nbImgGPU,'disableGPU',disableGPU);
    t=toc;fprintf(['\b\b\b (' num2str(t) 's) / ']);
    
%% FDL construction
    tic;fprintf(['FDL construction...']);
     if(useGPU)
         F = ComputeFDL_SuperRes_gpu(F, U, V, D, lambda, SRFact, dwnFilt ...
             ,'SRFactPrecond', SRFactPrecond....
             ,'deltaUV', deltaUV...
             ,'hexaSampling', hexaSampling...
             ,'lambdaColor', lambdaColor...
             ,'useAngPrecond', useAngPrecond...
             ,'sigAng', sigAng...
             ,'Px', Px, 'Py', Py...
             ,'ViewWeights', ViewWeights...
             ,'outputMode', outputModeCur...
             ,'freqBlockSz', freqBlockSz...
             );
    else
        error('CPU version of FDL super-resolution method not implemented yet');
    end
    t=toc;fprintf(['\b\b\b (' num2str(t) 's) / ']);
    
    if(it==numIterMax)
        fprintf('\n');
        break;
    end
    
%% Inverse Fourier transform and Square to Hexagonal resampling:
    tic;fprintf(['Inverse FFT...']);
    F = ifftImgs(F, [], 'dimsOrderIn', [4 3 1 2], 'dimsOrderOut', [1 2 3 4], 'hexaSampling',hexaSampling,'nbImgGPU',nbImgGPU,'disableGPU',disableGPU);
    t=toc;fprintf(['\b\b\b (' num2str(t) 's)\n']);
    
%% Replace known pixels
    if(computeMSE)
        IterInfo.MSE(it) = mean((ViewsRef(:)-real(F(:))).^2);
    end
    if(computeVideo && mod(it-1,seqIterStep)==0)
        IterInfo.seqIter(:,:,:,ceil(it/seqIterStep)) = F(:,:,:,ceil(numImages/2));
    end
    
    if(it<numIterMax)
        if(nonNegative)
            F = bsxfun(@times, 1-Mask, max(real(F),0)) + Views;
        else
            F = bsxfun(@times, 1-Mask, F) + Views;
        end
    end
    
end