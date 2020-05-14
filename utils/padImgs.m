% Performs padding (and windowing of the padded borders) for a set of images.
%
%-Inputs:
%   - Imgs:  Set of images given as a 4D array.
%            By default the dimensions must be in the order: 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color components, 4.Images,
%            but any permutation of the dimensions can be used by specifying the 'dimsOrderOut' optional argument.
%   - borderParams (optional): object of class BorderParams defining the parameters for padding and windowing of padded borders. The properties used are:
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
%       - disableGPU: Set to true to disable GPU, even if a compatible GPU is available (default=false).
%       - nbImgGPU: Number of Images treated in parallel for GPU computation (default=32, may need adjustements depending on gpu memory available).
%
%-Output:
%  - ImgsPad: Array of padded and transformed images.
%             By default the dimensions are in the order: 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color components, 4.Images,
%             but any permutation of the dimensions can be obtained by specifying the 'dimsOrderOut' optional argument.

% See also : fftImgs, BorderParams

function ImgsPad = padImgs(Imgs, borderParams, varargin)
p = inputParser;
addParameter (p,'dimsOrderIn',      [1 2 3 4],  @(x)isnumeric(x)&&numel(x)==4&&isequal(sort(x(:)),[1;2;3;4]));
addParameter (p,'dimsOrderOut',     [],  		@(x)isnumeric(x)&&numel(x)==4&&isequal(sort(x(:)),[1;2;3;4]));
addParameter (p,'disableGPU',   false,  @islogical);
addParameter (p,'nbImgGPU',     32,     @(x)isnumeric(x)&&x>0);
parse(p,varargin{:});
dimsOrderIn = p.Results.dimsOrderIn(:);
if(isempty(p.Results.dimsOrderOut))
	dimsOrderOut = dimsOrderIn;
else
	dimsOrderOut = p.Results.dimsOrderOut(:);
end
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

useWindow = ~strcmp(borderParams.window,'none');
useGPU = parallel.gpu.GPUDevice.isAvailable && ~disableGPU;

%% Transform input images (include padding and windowing).
%Prepare image size and windowing data.
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

%% Apply padding/windowing
ImgsPad = zeros([fullResY,fullResX,nChan,numImgs],class(Imgs));
if(useGPU)
    for l=1:nbImgGPU:numImgs
        idImgs = l:min(l+nbImgGPU-1,numImgs);        
        Tmp = padarray(gpuArray(Imgs(:,:,:,idImgs)), [padSizeY padSizeX], borderParams.padVal);
        if(exceedPadding), Tmp = Tmp(cropIdsY,cropIdsX,:,:); end
        if(useWindow), Tmp = bsxfun(@times, Tmp, W);end
        ImgsPad(:,:,:,idImgs) = gather(Tmp);
    end
    clear Tmp
else
    for l=1:numImgs
        for ch=1:nChan
            Tmp = padarray(Imgs(:,:,ch,l), [padSizeY padSizeX],borderParams.padVal);
            if(exceedPadding), Tmp = Tmp(cropIdsY, cropIdsX);end
            if(useWindow), Tmp = bsxfun(@times, Tmp, W);end
            ImgsPad(:,:,ch,l) = Tmp;
        end
    end
    clear Tmp
end

ImgsPad = permute(ImgsPad,dimsOrderOut);

end
