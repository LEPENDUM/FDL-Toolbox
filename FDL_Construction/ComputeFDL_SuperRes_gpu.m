%Compute Fourier Disparity Layers with spatial super-resolution using GPU,
%and reconstruct either the super-resolved views, the filtered low
%resolution views (with super-resolution and downsampling) or the high
%resolution FDL.
%
%----------- INPUTS -----------
%  -ViewsFFT:
%       Fourier Transform of input views with dimensions: 1.Views, 2.Color channels, 3.Vertical axis (Y), 4.Horizontal axis (X).
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
%  -SRFact:
%       Super-resolution factor: SRFact(1) = horizontal factor.
%                                SRFact(2) = vertical factor (if not defined, same as horizontal factor).
%
%  -dwnFilt:
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
%       -freqBlockSz (default=1024):
%           Number of frequencies processed in parallel.
%           It may be adjusted to prevent 'out of memory' depending on available GPU memory, number of views, and number of layers.
%       -outputMode:
%           -'LR': return the low resolution light field resulting from super-resolution and downsampling (i.e. low resolution FDL filtering preserving spatial aliasing).
%           -'HR': return the reconstructed High resolution light field (no residual back-projection).
%           -'HR-Backproj': return the reconstructed High resolution light field with residual back-projection.
%           -'HR-FDL' (default): return the high resolution FDL model.
%
%----------- OUTPUT -----------
%  -OutputFT: Output in the Fourier domain as specified by the 'outputMode' argument.
%             For a FDL output ('HR-FDL' mode), the dimensions are in the order : 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color channels, 4.Layers.
%             For views output, the dimensions are in the same order as the input: 1.Views, 2.Color channels, 3.Vertical axis (Y), 4.Horizontal axis (X).

function OutputFT = ComputeFDL_SuperRes_gpu(ViewsFFT, U, V, D, lambda, SRFact, dwnFilt, varargin)

p = inputParser;
addParameter (p,'outputMode',       'HR-FDL',   @(x)any(validatestring(x,{'LR','HR','HR-Backproj','HR-FDL'})));
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

parse(p,varargin{:});
outputMode = p.Results.outputMode;
lambdaColor = p.Results.lambdaColor;
hexaSampling = p.Results.hexaSampling;
deltaUV = p.Results.deltaUV;
useAngPrecond = p.Results.useAngPrecond;
sigAng = p.Results.sigAngular;
ViewWeights = p.Results.ViewWeights;
freqBlockSz = p.Results.freqBlockSz;
SRFactPrecond = p.Results.SRFactPrecond;
Px = p.Results.Px;
Py = p.Results.Py;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DeltaU = deltaUV(1);
if(numel(deltaUV)==1)
    DeltaV = DeltaU;
else
    DeltaV = deltaUV(2);
end

ViewsFFT = single(ViewsFFT);

outputFDL = false;
switch outputMode
    case 'LR'
        HR_Output = false;
        doBackProj = false;
    case 'HR'
        HR_Output = true;
        doBackProj = false;
    case 'HR-Backproj'
        HR_Output = true;
        doBackProj = true;
    case 'HR-FDL'
        outputFDL = true;
        HR_Output = true;
        doBackProj = false;
end

%% Set scale factor parameters
SRFact = ceil(SRFact);
SRFactX = SRFact(1);
if(numel(SRFact)==1)
    SRFactY = SRFactX;
elseif(numel(SRFact)==2)
    SRFactY = SRFact(2);
else
    error('Input argument SRFact must have 1 or 2 elements');
end
%Set preconditioning scale (may be different from Super-resolution scale e.g. to disable spatial preconditioning / to use preconditioning on already upsampled input).
SRFactPrecond = ceil(SRFactPrecond);
if(isempty(SRFactPrecond))
    SRFactPrecondX = SRFactX;
    SRFactPrecondY = SRFactY;
else
    SRFactPrecondX = SRFactPrecond(1);
    if(numel(SRFactPrecond)>1)
        SRFactPrecondY = SRFactPrecond(2);
    else
        SRFactPrecondY = SRFactPrecondX;
    end
end
useSpatPrecond = ( SRFactPrecondX>1 || SRFactPrecondY>1 );

%% Set downsampling filter parameters
if(iscell(dwnFilt))
    switch dwnFilt{1}
        case 'exactGauss'
            dwnFiltType=0;
            %gaussian blur parameter (relative to the input low resolution image size).
            sigLR = dwnFilt{2};
        case 'box'
            dwnFiltType=1;
            dwnFiltKernel = ones(SRFactY,SRFactX);
            dwnFiltKernelOX=1;dwnFiltKernelOY=1;
        case 'custom'
            dwnFiltType=1;
            %array representing the downsampling kernel (relative to the target high resolution image size).
            dwnFiltKernel = dwnFilt{2};
            %index of the kernel's origin (at the center by default).
            if(numel(dwnFilt)>2)
                dwnFiltKernelOX = dwnFilt{3};
                dwnFiltKernelOY = dwnFilt{4};
            else
                dwnFiltKernelOX = ceil((size(dwnFiltKernel,2)+1)/2);
                dwnFiltKernelOY = ceil((size(dwnFiltKernel,1)+1)/2);
            end
        otherwise
            error('unknown downsampling filter type');
    end
else
    dwnFiltType=0;
    sigLR=dwnFilt;
end

%% Set dimension parameters
numDisp = size(D,2);
numInitViews = size(ViewsFFT,1);
nChan = size(ViewsFFT,2);
numFreqX = size(ViewsFFT,4);
numFreqY = size(ViewsFFT,3);
numCoeff = SRFactX*SRFactY;
numDispSR = numCoeff*numDisp;
if(HR_Output)
    even_fft = 1-mod([numFreqY numFreqX].*[SRFactY SRFactX],2);
else
    even_fft = 1-mod([numFreqY numFreqX],2);
end
nChanParams = size(D,3);
if(nChanParams == 1)
    nChanParams = size(U,3);
end
if(nChanParams~=1 && nChanParams~=nChan)
    error('The number of channels for the calibration parameters should be either 1 or equal to the number of channels of the input data.');
end
if( ~isequal(size(U),size(V)) || (nChanParams~=1 && size(U,3)~=nChanParams) || (nChanParams~=1 && size(V,3)~=nChanParams) )
    error('Number of channel mismatch for U, V and D parameters');
end

if(doBackProj)
    %For residual back-projection, the residual deblurring can't be
    %disparity-dependant => If disparity-depenendant blur (i.e. angular blur)
    %is used, we choose the lowest disparity for residual deblurring (smallest deblurring).
    if(sigAng>0)
       [~,idMinD_disp]=min(sum(D.^2,3));
    else
        idMinD_disp=1;
    end
end


%% Compute frequency regions for Fourier domain convolution (equivalent to spatial domain multiplication by mask of unknown pixels).
%if(HR_Output)
    [wx,wy,xC,yC,xRanges,yRanges,coeffs,centralFreqsId]=GenerateFrequencies_SR(numFreqX,numFreqY,[SRFactX SRFactY],hexaSampling);
    %Coefficients of the Fourier coefficients of the Super-resolution mask (all equal to 1 for square resampling or hexaSampling=2)
    if(outputFDL)
        coeffs = reshape(gpuArray(coeffs),1,numCoeff,1);
    else
        coeffs = reshape(gpuArray(coeffs),1,1,1,numCoeff,1);
    end
%else
%    [wx,wy,xC,yC,~,~,~,centralFreqsId]=GenerateFrequencies_SR(numFreqX,numFreqY,[SRFactX SRFactY],hexaSampling);
%end
wx = permute(gpuArray(reshape(wx,[],numCoeff)),[5 4 3 2 1]);
wy = permute(gpuArray(reshape(wy,[],numCoeff)),[5 4 3 2 1]);
numFreqProc = numFreqY*(xC-1)+yC; %process half frequencies (remaining freqencies obtained by symmetry).


%% Prepare filter data:
%spatial downsampling filter
if(dwnFiltType==0)
    sig=sigLR(:).*SRFact(:);%gaussian blur std.dev. proportional to the image scale.
    sigFactX = 2*(pi*sig(1)).^2;
    if(numel(sig)==1)
        sigFactY = sigFactX;
    else
        sigFactY = 2*(pi*sig(2)).^2;
    end
    useDwnFilt = any(sig>0);
else
    [shDwnX, shDwnY]=meshgrid(1-dwnFiltKernelOX:size(dwnFiltKernel,2)-dwnFiltKernelOX, dwnFiltKernelOY-1:-1:dwnFiltKernelOY-size(dwnFiltKernel,1));
    shDwnX=shDwnX(:);
    shDwnY=shDwnY(:);
    cDwn=dwnFiltKernel(:)/sum(dwnFiltKernel(:));
    useDwnFilt = sum(abs(cDwn))>0;
end
%angular blur filter.
sigFactAng = 2*(pi*sigAng).^2;

%spatial preconditioning filter.
if(useSpatPrecond)
    [shPreX, shPreY]=meshgrid(0:SRFactPrecondX-1,0:-1:1-SRFactPrecondY);
    shPreX=shPreX(:);
    shPreY=shPreY(:);
    cPre=ones(numCoeff,1)/numCoeff;
end



%% Set FDL model parameters
U = single(reshape(U,numInitViews,1,nChanParams))*SRFactX;
V = single(reshape(V,numInitViews,1,nChanParams))*SRFactY;
D = single(reshape(D,1,numDisp,nChanParams));

if(~exist('Px','var') || ~exist('Py','var') || isempty(Px) || isempty(Py))
    Px = pagefun( @mtimes, gpuArray(single(U)), gpuArray(single(D)) );
    Py = pagefun( @mtimes, gpuArray(single(V)), gpuArray(single(D)) );
else
    nChanParams = size(Px,3);
    if(~isequal(size(Px),size(Py)) || size(Px,1)~=numInitViews || (nChanParams~=1 && nChanParams~=nChan) )
        error('Dimension mismatch for parameter matrices Px and Py.');
    end
    Px = gpuArray(single(Px))*SRFactX;
    Py = gpuArray(single(Py))*SRFactY;
end
%Take the same view positions and layer disparity parameters for all color components for final reconstruction (->reduce chromatic aberrations).
if(HR_Output)
    DReco = mean(D,3);
    PxReco = mean(U,3)*DReco;
    PyReco = mean(V,3)*DReco;
end

%View Weights parameters (weight the different views in the FDL construction).
useViewWeights = ~isempty(ViewWeights);
if(useViewWeights)
    if(numel(ViewWeights)~=numInitViews)
        error('If ''ViewWeights'' is not empty, the number of elements should be equal to the number of views.');
    end
    ViewWeights = reshape(gpuArray(single(ViewWeights)),1,[]).^2;
end


DeltaU = DeltaU*SRFactX;%adjust horizontal baseline to target resolution (original parameters ).
DeltaV = DeltaV*SRFactY;

%% Adjustment of the regularization parameter to keep the results consistent with repect to number of input views and number of layers.
lambdaL2 = 5e-3*lambda*numInitViews*numDisp/numCoeff;
lambdaRGBReg = lambdaColor*numInitViews*numDisp*numCoeff/175; %regularization parameter for the 2nd Order view regularization.
%% Prepare matrices that do not depend on frequency.
%L2 Regularization matrix.
if(lambdaRGBReg>0)
    Reg = lambdaL2*gpuArray(eye(numDispSR*nChan)); 
else
    Reg = lambdaL2*gpuArray(eye(numDispSR)); 
end
dupliChan = gpuArray(reshape(eye(nChan),[1 nChan 1 nChan]));%Matrix used to apply kronecker product with nChan x nChan identity.
dupliCoeffs = gpuArray(reshape(eye(numCoeff),1,numCoeff,1,numCoeff));%Matrix used to apply kronecker product with numCoeff x numCoeff identity.
%Color regularization parameters:
% 1.Discrete spatial laplacian filter parameters.
gradShifts = [ 0 -1  1  0  0 ; -1  0  0  1  0];
gradCoeffs = [-1 -1 -1 -1  4];
% 2.Matrix for evaluating color difference constancy.
GTGBase = gpuArray(single([-2 .5 1.5;2 0 -2]));
GTGBase = reshape(gpuArray(single(GTGBase'*GTGBase)),1,3,1,3);


%% FDL Reconstruction
if(HR_Output)
    if(outputFDL)
        OutputFT = zeros(numDisp,numCoeff,nChan,numFreqX*numFreqY,'single');%output = High resolution FDL.
    else
        OutputFT = zeros(numInitViews,nChan,numCoeff,numFreqX*numFreqY,'single');%output = High resolution views.
    end
else
    OutputFT = zeros(size(ViewsFFT),'single');%output = Low resolution views (after super-resolution followed by downsampling).
end

for freqSt=1:freqBlockSz:numFreqProc
    freqEnd = min(freqSt+freqBlockSz-1,numFreqProc);
    freqList = freqSt:freqEnd;
    
    ViewsFFT_gpu = gpuArray(ViewsFFT(:,:,freqList));
    
    %Blur Model filters
    if(useDwnFilt)%Spatial blur model
        if(dwnFiltType==0)%Exact gaussian blur model
            FiltModel = exp(-(sigFactX*wx(:,:,:,:,freqList).^2+sigFactY*wy(:,:,:,:,freqList).^2));
        else              %Discrete blur kernel
            FiltModel = sum(bsxfun(@times,exp(2i*pi*(bsxfun(@times,shDwnX,wx(:,:,:,:,freqList))+bsxfun(@times,shDwnY,wy(:,:,:,:,freqList)))),cDwn),1);
        end
    else
        FiltModel = ones(1,1,1,numCoeff,1);
    end
    if(sigAng>0)%Angular blur model (depth-dependent gaussian)
        FiltModel = bsxfun(@times,FiltModel,...
            exp(bsxfun(@times,D.^2,-(sigFactAng*DeltaU^2*wx(:,:,:,:,freqList).^2+sigFactAng*DeltaV^2*wy(:,:,:,:,freqList).^2))));
    end
    
    %Spatial preconditionning filter.
    if(useSpatPrecond)
        %Default Spatial preconditioning (without correction for blur model) ==> magnitude of nearest interpolation filter (=square root of bilinear filter).
        FiltPrecond = abs(sum(bsxfun(@times,exp(2i*pi*(bsxfun(@times,shPreX,wx(:,:,:,:,freqList))+bsxfun(@times,shPreY,wy(:,:,:,:,freqList)))),cPre),1));
        %Spatial Super-Resolution Preconditioning correction accounting for spatio-angular blur model.
        if(useDwnFilt || sigAng>0)
            FiltModelSq = abs(FiltModel).^2;
        %3.Optimised to apply upsampling followed by deconvolution in the worst case (i.e. single view).
            FiltPrecond = bsxfun(@times,FiltPrecond,sqrt(lambdaL2./(bsxfun(@times,lambdaL2+FiltModelSq,1-sum(bsxfun(@times,FiltPrecond.^2,FiltModelSq./(lambdaL2+FiltModelSq)),4)))));
            clear FiltModelSq
        end
    else
        FiltPrecond=1;
    end
    
    %Angular preconditionning filter.
    if(useAngPrecond)
        FiltPrecond = bsxfun(@times,FiltPrecond, sinc(DeltaU*bsxfun(@times,D,wx(:,:,:,:,freqList))).*sinc(DeltaV*bsxfun(@times,D,wy(:,:,:,:,freqList))));
    end
    
    %Combine preconditionning and blur model.
    A = FiltPrecond.*FiltModel;
    if(~doBackProj)
        clear FiltModel
    end
    if(~HR_Output)
        clear FiltPrecond
    end
    
    %Construct FDL model matrix.
    A = bsxfun(@times,A, exp( 2i*pi*(bsxfun(@times,Px,wx(:,:,:,:,freqList)) + bsxfun(@times,Py,wy(:,:,:,:,freqList)))));
    %Concatenated matrices for related frequencies in the convolution (Fourier domain convolution = spatial pixel grid masking).
    A = reshape(permute(A,[1 2 4 3 5]), numInitViews,numDispSR,nChanParams,[]);
    
    AT = pagefun(@ctranspose,A);
    if(useViewWeights)
        AT = bsxfun(@times,AT,ViewWeights);
    end
    ATA = pagefun(@mtimes,AT,A);
    
    
    if(HR_Output && ~doBackProj)
        clear A
    end
    
    %Premulitply input views data with AT to free AT memory before caling linear least square solver.
    ViewsFFT_gpu = pagefun(@mtimes,AT, reshape(ViewsFFT_gpu,numInitViews,1,nChan,[]) );
    clear AT
    
    
    if(lambdaRGBReg>0)
        %Build color regularization matrix:
        %1.Process discrete laplacian filter.
        Grad2 = abs(sum(bsxfun(@times,gradCoeffs,exp( bsxfun(@times,2i*pi*wx(:,:,:,:,freqList),gradShifts(1,:)) + bsxfun(@times,2i*pi*wy(:,:,:,:,freqList),gradShifts(2,:)) )),2)).^2;
        Grad2 = bsxfun(@times, gpuArray(single(eye(numDisp))), Grad2);
        %2.Form differences of the filtered color components.
        Grad2 = reshape(bsxfun(@times, dupliCoeffs, reshape(Grad2,numDisp,1,numDisp,numCoeff,[])), numDispSR, numDispSR,[]);
        Grad2 = reshape(bsxfun(@times, GTGBase, reshape(Grad2,numDispSR,1,numDispSR,1,[])), nChan*numDispSR, nChan*numDispSR,[]);
        RegFull = Reg + lambdaRGBReg*Grad2;
        
        %Convert ATA to block diagonal matrix and vectorize ViewsFFT_gpu along color components to make matrix sizes compatible with color regularization matrix.
        ATA = reshape(bsxfun(@times, reshape(ATA,numDispSR,1,numDispSR,nChanParams,[]), dupliChan), numDispSR*nChan, numDispSR*nChan, []);
        ViewsFFT_gpu = reshape(ViewsFFT_gpu,numDispSR*nChan,1,[]);

        %Solve Ridge Regression (linear least squares with Tikhonov Regularization)
        FDL_gpu = ...
            reshape(...
            pagefun(@mldivide, bsxfun(@plus, ATA, RegFull), ViewsFFT_gpu),...
            numDispSR,1,nChan,[]);
        
    else
        %Solve Ridge Regression (linear least squares with Tikhonov Regularization)
        FDL_gpu = pagefun(@mldivide, bsxfun(@plus, ATA, Reg), ViewsFFT_gpu);
    end
    clear ATA
    
    if(HR_Output && ~outputFDL)
%Output = high resolution views (with or without back-projection).
        if(doBackProj)
            R = reshape(gpuArray(ViewsFFT(:,:,freqList)),numInitViews,1,nChan,1,[]) - reshape(pagefun(@mtimes,A,FDL_gpu),numInitViews,1,nChan,1,[]);
        end
        
        A = exp( 2i*pi*(bsxfun(@times,PxReco,wx(:,:,:,:,freqList)) + bsxfun(@times,PyReco,wy(:,:,:,:,freqList))) );
        A = bsxfun(@times, FiltPrecond, A);
        
        FDL_gpu = permute(reshape(FDL_gpu,numDisp,numCoeff,nChan,[]),[1 5 3 2 4]);
        FDL_gpu = numCoeff*bsxfun(@times,FDL_gpu,coeffs);%Adjustments necessary for the High resolution output (factor numCoeff / apply coeffs of super-resolution mask).
        
        if(doBackProj)
            RecHR = pagefun(@mtimes,A,FDL_gpu);
            RecHR(:,:,:,centralFreqsId,:) = RecHR(:,:,:,centralFreqsId,:) + bsxfun(@times,R, numCoeff./FiltModel(:,idMinD_disp,:,centralFreqsId,:));%with deconvolution of the residual (only spatial) after upsampling.
            OutputFT(:,:,:,freqList) = gather(RecHR);
        else
            OutputFT(:,:,:,freqList) = gather(pagefun(@mtimes,A,FDL_gpu));
        end
    elseif(outputFDL)
%Output = high resolution FDL.
        FDL_gpu = reshape(FDL_gpu,numDisp,numCoeff,nChan,[]);
        FDL_gpu = numCoeff*bsxfun(@times,FDL_gpu,coeffs);%Adjustments necessary for the High resolution output (factor numCoeff / apply coeffs of super-resolution mask).
        FDL_gpu = bsxfun(@times,FDL_gpu,permute(FiltPrecond,[2,4,3,5,1]));
        OutputFT(:,:,:,freqList) = gather(FDL_gpu);
    else
%Output = low resolution views.
        OutputFT(:,:,freqList) = gather(pagefun(@mtimes,A,FDL_gpu));
    end
end


clear A ViewsFFT_gpu FDL_gpu


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Rearrange frequencies into the 2D grid                     %
%              + Find other frequencies with symmetry                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(HR_Output && numCoeff>1)
%For final reconstruction, the output is in high resolution => the frequencies must be arranged in the frequency spectrum for high resolution.
    ids = zeros(numCoeff,numFreqProc);
    ids_sym = zeros(numCoeff,numFreqProc);
    for i=1:numCoeff
        [idX,idY]=meshgrid(xRanges(:,i),yRanges(:,i));
        ids(i,:)=sub2ind([numFreqY,numFreqX].*[SRFactY SRFactX],idY(1:numFreqProc),idX(1:numFreqProc));

        idX_sym = mod(SRFactX*numFreqX - idX+even_fft(2), SRFactX*numFreqX)+1;
        idY_sym = mod(SRFactY*numFreqY - idY+even_fft(1), SRFactY*numFreqY)+1;
        ids_sym(i,:)=sub2ind([numFreqY,numFreqX].*[SRFactY SRFactX],idY_sym(1:numFreqProc),idX_sym(1:numFreqProc));
    end
    ids=reshape(ids,1,[]);
    ids_sym=reshape(ids_sym,1,[]);
    
    if(outputFDL)
        OutputFT = reshape(permute(OutputFT,[2 4 3 1]),numCoeff*numFreqX*numFreqY,nChan,numDisp);
        OutputFT(ids,:,:) = OutputFT(1:numFreqProc*numCoeff,:,:);
        OutputFT(ids_sym,:,:) = conj(OutputFT(ids,:,:));
        OutputFT = reshape(OutputFT,SRFactY*numFreqY,SRFactX*numFreqX,nChan,numDisp);
    else
        OutputFT = reshape(OutputFT,numInitViews,nChan,numCoeff*numFreqX*numFreqY);
        OutputFT(:,:,ids) = OutputFT(:,:,1:numFreqProc*numCoeff);
        OutputFT(:,:,ids_sym) = conj(OutputFT(:,:,ids));
        OutputFT = reshape(OutputFT,numInitViews,nChan,SRFactY*numFreqY,SRFactX*numFreqX);
    end
    
else
%For an output without super-resolution (numCoeff==1) or with super-resolution followed by downsampling
%=> the frequency spectrum is arranged the same way as input.
    if(outputFDL)
        OutputFT = reshape(OutputFT,numDisp,nChan,numFreqY,numFreqX);
    else
        OutputFT = reshape(OutputFT,numInitViews,nChan,numFreqY,numFreqX);
    end
    
    OutputFT(:,:,yC+1:end,xC:numFreqX) = conj(OutputFT(:,:,yC-1:-1:1+even_fft(1),xC:-1:1+even_fft(2)));
    OutputFT(:,:,1+even_fft(1):yC,xC+1:numFreqX) = conj(OutputFT(:,:,end:-1:yC,xC-1:-1:1+even_fft(2)));
    if(even_fft(1))
        OutputFT(:,:,1,xC+1:numFreqX) = conj(OutputFT(:,:,1,xC-1:-1:1+even_fft(2)));
    end
    
    if(outputFDL)
        OutputFT = permute(OutputFT,[3 4 2 1]);
    end
end

end