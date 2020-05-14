% Performs 2D Fast Fourier transform of a set of images.
% Padding parameters can be specified to perform a padding (and windowing
% of the padded borders) before applying the Fourier transform.
%
%-Inputs:
%   - Imgs: Set of images in pixel domain given as a 4D array.
%           By default the dimensions must be in the order: 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color components, 4.Images,
%           but any permutation of the dimensions can be used by specifying the 'dimsOrderOut' optional argument.
%   - borderParams (optional): Object of class BorderParams defining the parameters for padding and windowing of padded borders. The properties used are:
%       - L: padding on the left side (default=0).
%       - R: padding on the right side (default=0).
%       - T: padding on the top (default=0).
%       - B: padding on the bottom (default=0).
%       - padVal: padding value: either numeric value or 'replicate', 'symmetric' or 'circular' (default='symmetric').
%       - window:
%           - 'none' -> no windowing of padded borders.
%           - 'hann' (default) -> use hann window type for padded borders.
%           - 'linear' -> use linear window type for padded borders.
%  -varargin: Additional optional arguments given as (Name,Value) pairs:
%       - dimsOrderIn (default=[1 2 3 4]): Order of the input dimensions indicated as a permutation of the vector [1 2 3 4] where:
%           1 represents the vertical axis (Y).
%           2 represents the horizontal axis (X).
%           3 represents the dimension of color components.
%           4 represents the dimension of the image set.
%         e.g. use dimsOrderIn=[4 3 1 2] for an input array of size (#Images x #Color Components x Vertical Resolution x Horizontal Resolution).
%       - dimsOrderOut (equal to dimsOrderIn by default): Order of the output dimensions (same format as 'dimsOrderIn').
%       - hexaSampling : Value indicating the sampling format of the input images:
%           -0 (default): Square sampling (use normal 2D FFT).
%           -1: Hexagonal sampling where rows with odd indices are shifted to the left by half a pixel (where indices start at 1 for the top row of the images excluding padded borders).
%           -2: Hexagonal sampling where rows with even indices are shifted to the left by half a pixel (where indices start at 1 for the top row of the images excluding padded borders).
%           Note: For hexagonal samplings with shifting of the columns instead of rows, use the 'dimsOrderIn' and 'dimsOrderOut' arguments to permute horizontal and vertical dimensions.
%       - disableGPU: Set to true to disable GPU, even if a compatible GPU is available (default=false).
%       - nbImgGPU: Number of Images treated in parallel for GPU computation (default=32, may need adjustements depending on gpu memory available).
%
%-Output:
%  - ImgsFFT: Array of padded and transformed images.
%             By default the dimensions are in the order: 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color components, 4.Images,
%             but any permutation of the dimensions can be obtained by specifying the 'dimsOrderOut' optional argument.

% See also : ifftImgs, padImgs, BorderParams

function ImgsFFT = fftImgs(Imgs, borderParams, varargin)
p = inputParser;
addParameter (p,'dimsOrderIn',      [1 2 3 4],  @(x)isnumeric(x)&&numel(x)==4&&isequal(sort(x(:)),[1;2;3;4]));
addParameter (p,'dimsOrderOut',     [],  		@(x)isnumeric(x)&&numel(x)==4&&isequal(sort(x(:)),[1;2;3;4]));
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
hexaSampling = p.Results.hexaSampling;
disableGPU = p.Results.disableGPU;
nbImgGPU = p.Results.nbImgGPU;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('borderParams','var') || isempty(borderParams)), borderParams = BorderParams();end
if(~isa(borderParams,'BorderParams'))
    error('Input argument ''borderParams'' should be either empty or a BorderParams object.');
end

padSizeX = max(borderParams.L, borderParams.R);
padSizeY = max(borderParams.T, borderParams.B);
exceedPadding = borderParams.L~=borderParams.R || borderParams.T~=borderParams.B;
if(exceedPadding)
    cropL = max(0,borderParams.R-borderParams.L);
    cropR = max(0,borderParams.L-borderParams.R);
    cropT = max(0,borderParams.B-borderParams.T);
    cropB = max(0,borderParams.T-borderParams.B);
end

if(hexaSampling>0)
    %Select first row index to shift so that hexaSampling defines the parity of row indices
    %where indices start at 1 for the top row of the images excluding padded borders.
    shiftRowStart = 2 - mod(hexaSampling + mod(borderParams.T,2), 2);%swap values 1 and 2 if borderParams.T is odd.
end

useWindow = ~strcmp(borderParams.window,'none') && max(padSizeX,padSizeY)>0;
useGPU = parallel.gpu.GPUDevice.isAvailable && ~disableGPU;

%% Transform input images (include padding and windowing).
%Prepare dimensions variables.
dimsInPos(dimsOrderIn)=[1 2 3 4];
resY = size(Imgs,dimsInPos(1));
resX = size(Imgs,dimsInPos(2));
nChan = size(Imgs,dimsInPos(3));
numImgs = size(Imgs,dimsInPos(4));
fullResX = resX + borderParams.L + borderParams.R;
fullResY = resY + borderParams.T + borderParams.B;

%Prepare window for padded borders
if(useWindow)
    W = single(GenerateWindowIntensity(resX+2*padSizeX,resY+2*padSizeY,padSizeX,padSizeY,borderParams.window));
    if(useGPU), W = gpuArray(W);end
    if(exceedPadding), W = W(1+cropT:end-cropB, 1+cropL:end-cropR);end
end

%Precompute index sets.
if(exceedPadding)
	%Indices for cropping (crop of 1 px at most in the case of assymetric padding).
	cropIdsY = 1+cropT:resY+2*padSizeY-cropB;
	cropIdsX = 1+cropL:resX+2*padSizeX-cropR;
end
if(hexaSampling~=0)
	%odd or even indices selection for hexagonal sampling
	hexIds = shiftRowStart:2:fullResY;
end

%Precompute Fourier domain shift operation for hexagonal sampling.
if(hexaSampling>0)
    xC = ceil((fullResX+1)/2);
    wx = single(([1:fullResX]-xC)/(fullResX));
    hexaShift = reshape(exp(-1i*pi*wx(1,:)),1,[]);%precompute right shift by 1/2 pixel (used to compensate for the left shift of input rows in hexagonal sampling).
    if(useGPU), hexaShift = gpuArray(hexaShift);end
    clear wx xC
end

%% Apply transform.
ImgsFFT = zeros(fullResY,fullResX,nChan,numImgs,'single');
Imgs = permute(Imgs,dimsInPos);
if(useGPU)
    for l=1:nbImgGPU:numImgs
        idImgs = l:min(l+nbImgGPU-1,numImgs);%Indices of selected images to process in parallel
        %Padding
        Tmp = padarray(gpuArray(Imgs(:,:,:,idImgs)), [padSizeY padSizeX], borderParams.padVal);
        if(exceedPadding), Tmp = Tmp(cropIdsY,cropIdsX,:,:); end
        if(useWindow), Tmp = bsxfun(@times, Tmp, W);end
        %Fourier Transform
        if(hexaSampling==0)
            %Square sampling -> normal 2D FFT.
            Tmp = fftshift(fftshift(fft(fft(Tmp,[],1),[],2),1),2);
        else
            %FFT assuming input data given with hexagonal sampling.
            Tmp = fftshift(fft(Tmp,[],2),2); %FFT of rows.
			Tmp(hexIds,:,:,:) = bsxfun(@times, Tmp(hexIds,:,:,:), hexaShift); %Shift odd or even rows by half a pixel.
            Tmp = fftshift(fft(Tmp,[],1),1); %FFT of columns.
        end
        ImgsFFT(:,:,:,idImgs) = gather(Tmp);
    end
    
    clear Tmp
else
    for l=1:numImgs
        for ch=1:nChan
            %Padding
            Tmp = padarray(Imgs(:,:,ch,l), [padSizeY padSizeX],borderParams.padVal);
            if(exceedPadding), Tmp = Tmp(cropIdsY, cropIdsX);end
            if(useWindow), Tmp = bsxfun(@times, Tmp, W);end
            %Fourier Transform
            if(hexaSampling==0)
                %Square sampling -> normal 2D FFT
                ImgsFFT(:,:,ch,l) = fftshift(fft2(Tmp));
            else
                %FFT assuming input data given with hexagonal sampling.
                Tmp = fftshift(fft(Tmp,[],2),2);%FFT of rows.
                Tmp(hexIds,:) = bsxfun(@times, Tmp(hexIds,:), hexaShift); %Shift odd or even rows by half a pixel.
                ImgsFFT(:,:,ch,l) = fftshift(fft(Tmp,[],1),1); %FFT of columns.
            end
        end
    end
    clear Tmp
end

ImgsFFT = permute(ImgsFFT, dimsOrderOut);

end
