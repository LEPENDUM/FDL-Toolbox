% Performs 2D Fast Inverse Fourier transform of a set of images.
% The function assumes that the signal to be recovered is real, in order to
% accelerate the computations using the symmetry of the Fourier transform
% in the real case.
% Border parameters can be specified to perform a crop in the spatial
% domain (may be needed if the signal was padded before forward fourier
% transform).
%
%-Required Input:
%   - Imgs: Set of images in the complex Fourier domain given as a 4D array.
%           By default the dimensions must be in the order: 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color components, 4.Images,
%           but any permutation of the dimensions can be used by specifying the 'dimsOrderOut' optional argument.
%           The number of elements of the horizontal dimension can be ceil((ResX+1)/2) instead ResX since horizontal symmetry is exploited
%           (assuming the output signal is real). In this case, the function assumes that 'ImgsFFT' contains the left half of the spectrum
%           and the full horizontal resolution must be specified as the 'fullResX' argument.
%   - borderParams (ptional): Object of class BorderParams defining the parameters for cropping the borders (previously padded before performing forward fft). The properties used are:
%       - L: crop on the left side (default=0).
%       - R: crop on the right side (default=0).
%       - T: crop on the top (default=0).
%       - B: crop on the bottom (default=0).
%  -varargin: Additional optional arguments given as (Name,Value) pairs:
%       - dimsOrderIn (default=[1 2 3 4]): Order of the input dimensions indicated as a permutation of the vector [1 2 3 4] where:
%           1 represents the vertical axis (Y).
%           2 represents the horizontal axis (X).
%           3 represents the dimension of color components.
%           4 represents the dimension of the image set.
%         e.g. use dimsOrderIn=[4 3 1 2] for an input array of size (#Images x #Color Components x Vertical Resolution x Horizontal Resolution).
%       - dimsOrderOut (equal to dimsOrderIn by default): Order of the output dimensions (same format as 'dimsOrderIn').
%       - fullResX: Full horizontal resolution (default=[], i.e. set automatically using size of 'ImgsFFT').
%                   'fullResX' must be specified if the input images 'ImgsFFT' only contains the left half of the spectrum.
%       - hexaSampling : Value indicating the sampling format for the output inverse transformed images:
%           -0 (default): Square sampling (use normal 2D inverse FFT).
%           -1: Hexagonal sampling where rows with odd indices are shifted to the left by half a pixel (where indices start at 1 for the top row of the images excluding padded borders).
%           -2: Hexagonal sampling where rows with even indices are shifted to the left by half a pixel (where indices start at 1 for the top row of the images excluding padded borders).
%           Note: For hexagonal samplings with shifting of the columns instead of rows, use the 'dimsOrderIn' and 'dimsOrderOut' arguments to permute horizontal and vertical dimensions.
%       - disableGPU: Set to true to disable GPU, even if a compatible GPU is available (default=false).
%       - nbImgGPU: Number of Images treated in parallel for GPU computation (default=32, may need adjustements depending on gpu memory available).
%
%-Output:
%  - Imgs: Array of Inverse transformed images.
%          By default the dimensions are in the order: 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color components, 4.Images,
%          but any permutation of the dimensions can be obtained by specifying the 'dimsOrderOut' optional argument.

% See also : fftImgs, BorderParams

function Imgs = ifftImgs(ImgsFFT, borderParams, varargin)
p = inputParser;
addParameter (p,'dimsOrderIn',      [1 2 3 4],  @(x)isnumeric(x)&&numel(x)==4&&isequal(sort(x(:)),[1;2;3;4]));
addParameter (p,'dimsOrderOut',     [],  		@(x)isnumeric(x)&&numel(x)==4&&isequal(sort(x(:)),[1;2;3;4]));
addParameter (p,'fullResX',         [],         @isnumeric);
addParameter (p,'hexaSampling',     0,          @(x)isnumeric(x)&&(x==0||x==1||x==2));
addParameter (p,'disableGPU',       false,      @islogical);
addParameter (p,'nbImgGPU',         32,         @(x)isnumeric(x)&&x>0);
parse(p,varargin{:});
dimsOrderIn = p.Results.dimsOrderIn(:);
if(isempty(p.Results.dimsOrderOut))
	dimsOrderOut = dimsOrderIn;
else
	dimsOrderOut = p.Results.dimsOrderOut(:);
end
fullResX = p.Results.fullResX;
hexaSampling = p.Results.hexaSampling;
disableGPU = p.Results.disableGPU;
nbImgGPU = p.Results.nbImgGPU;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('borderParams','var') || isempty(borderParams)), borderParams = BorderParams();end
if(~isa(borderParams,'BorderParams'))
    error('Input argument ''borderParams'' should be either empty or a BorderParams object.');
end

useGPU = parallel.gpu.GPUDevice.isAvailable && ~disableGPU;

if(hexaSampling>0)
    %Select first row index to shift so that hexaSampling defines the parity of row indices
    %where indices start at 1 for the top row of the images excluding padded borders.
    shiftRowStart = 2 - mod(hexaSampling + mod(borderParams.T,2), 2);%swap values 1 and 2 if borderParams.T is odd.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimsInPos(dimsOrderIn)=[1 2 3 4];
resY = size(ImgsFFT,dimsInPos(1));
if(isempty(fullResX))
    resX = size(ImgsFFT,dimsInPos(2));
else
    resX = fullResX;
end
xC = ceil((resX+1)/2);
nChan = size(ImgsFFT,dimsInPos(3));
numImgs = size(ImgsFFT,dimsInPos(4));

%Precompute index sets.
xIds = 1:xC;
if(hexaSampling~=0)
	%odd or even indices selection for hexagonal sampling
	hexIds = shiftRowStart:2:resY;
end

%Prepare shift operation for hexagonal sampling.
if(hexaSampling>0)
    wx = single(([xC:-1:1]-xC)/(resX));
    hexaShift = reshape(exp(-1i*pi*wx(1,:)),1,[]);%shift to the right by 1/2 pixel.
    if(useGPU), hexaShift = gpuArray(hexaShift);end
    clear wx
end

%% Apply inverse transform.
Imgs = zeros(resY,resX,nChan,numImgs,'single');
dimsInPos(dimsOrderIn)=[1 2 3 4];
ImgsFFT = permute(ImgsFFT,dimsInPos);
if(useGPU)
    for l=1:nbImgGPU:numImgs
        idImgs = l:min(l+nbImgGPU-1,numImgs);%Indices of selected images to process in parallel
        %Inverse Fourier Transform
        Tmp = gpuArray(ImgsFFT(:,xIds,:,idImgs));
        %Tmp = gpuArray(ImgsFFT(:,:,:,idImgs));
        if(hexaSampling==0)
            %Square sampling -> normal inverse 2D FFT.
            Tmp = ifft(conj(fliplr(ifft(ifftshift(Tmp,1),[],1))),resX,2,'symmetric');
        else
            %Inverse FFT with hexagonal output sampling.
            Tmp = conj( fliplr(ifft(ifftshift(Tmp,1),[],1)));
            Tmp(hexIds,:,:) = bsxfun(@times, Tmp(hexIds,:,:), hexaShift );
            Tmp = ifft(Tmp,resX,2,'symmetric');
        end
        Imgs(:,:,:,idImgs) = gather(Tmp);
    end
    clear Tmp
else
    for l=1:numImgs
        for ch=1:nChan
            %Inverse Fourier Transform
            if(hexaSampling==0)
                %Square sampling -> normal inverse 2D FFT.
                Imgs(:,:,ch,l) = ifft(conj(fliplr(ifft(ifftshift(ImgsFFT(:,xIds,ch,l),1),[],1))),resX,2,'symmetric');
            else
                %Inverse FFT with hexagonal output sampling.
                Tmp = conj(fliplr(ifft(ifftshift(ImgsFFT(:,xIds,ch,l),1),[],1)));
                Tmp(hexIds,:) = bsxfun(@times, Tmp(hexIds,:), hexaShift );
                Imgs(:,:,ch,l) = ifft(Tmp,resX,2,'symmetric');
            end
        end
    end
    clear Tmp
end

%Crop and permutation
Imgs = permute( Imgs(1+borderParams.T:end-borderParams.B, 1+borderParams.L:end-borderParams.R,:,:), dimsOrderOut);
