% Shift and Sum Refocusing for Light Fields using FFT.
% The method is similar to conventional shift and sum, but the shifting
% operation of each view is performed in the Fourier domain (faster when
% generating several refocused images, or when the input views are already
% given in the Fourier domain).
%
% -Required Inputs:
%   - Views : Set of views in a 4D array with dimensions in the order:
%     1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color components, 4.Images.
%   - U     : List of horizontal angular position of each view.
%   - V     : List of vertical angular position of each view.
%   - dRange  : disparity value for refocusing (if a vector is given, the
%   function generates a focal stack video containing the refocused images
%   with each refocus disparity in dRange).
%
% -Optional Inputs:
%   - padParams: object of class BorderParams defining the parameters for
%   padding and windowing of padded borders within the processing. The
%   padding is only performed for an input in the spatial domain (see
%   'inputSpatial' option). After refocusing, the padded borders are
%   cropped only if the 'outputSpatial' option is active.
%   The properties used from the BorderParams object are:
%       - L: padding/crop on the left side (default=0).
%       - R: padding/crop on the right side (default=0).
%       - T: padding/crop on the top (default=0).
%       - B: padding/crop on the bottom (default=0).
%       - padType: type of padding: 'replicate', 'symmetric' or 'circular' (default='symmetric').
%       - window:
%           - 'none' -> no windowing of padded borders.
%           - 'hann' (default) -> use hann window type for padded borders.
%           - 'linear' -> use linear window type for padded borders.
%  -varargin: Additional optional arguments given as (Name,Value) pairs:
%       - inputSpatial (default=true):
%           - true -> assumes input views are in spatial domain.
%           - false-> assumes input views are already in Fourier domain.
%       - outputSpatial (default=true):
%           - true -> return output images in spatial domain.
%           - false-> return output images in Fourier domain (only left half of the spectrum).
%       - fullResX: (only used if inputSpatial is false) Full horizontal resolution (default=[], i.e. set automatically using size of Views).
%                   'fullResX' must be specified if the input views ViewsFFT only contains the left half of the spectrum.
%       - hexaSampling : (only used if inputSpatial is true) Value indicating the sampling format of the input images:
%           -0 (default): Square sampling (use normal 2D FFT).
%           -1: Hexagonal sampling where rows with odd indices are shifted to the left by half a pixel (where indices start at 1 for the top row of the images excluding padded borders).
%           -2: Hexagonal sampling where rows with even indices are shifted to the left by half a pixel (where indices start at 1 for the top row of the images excluding padded borders).
%       - disableGPU: Set to true to disable GPU, even if a compatible GPU is available (default=false).
%       - nbImgGPU: Number of Images treated in parallel for GPU computation (default=32, may need adjustements depending on gpu memory available).
%       - nbRefocusGPU: Number of refocusing processing in parallel for GPU (default=10, may need adjustements depending on gpu memory available).
%
% -Output:
%   - I: 4D array of refocused images with dimensions in the order:
%     1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color components, 4.Images.

% See also : fftImgs, ifftImgs, BorderParams

function I = refocusFFT(Views, U, V, dRange, padParams, varargin)
p = inputParser;
addParameter (p,'inputSpatial',     true,       @islogical);
addParameter (p,'outputSpatial',	true,       @islogical);
addParameter (p,'fullResX',         [],         @isnumeric);
addParameter (p,'hexaSampling',     0,          @(x)isnumeric(x)&&(x==0||x==1||x==2));
addParameter (p,'disableGPU',       false,      @islogical);
addParameter (p,'nbImgGPU',         32,         @(x)isnumeric(x)&&x>0);
addParameter (p,'nbRefocusGPU',     10,         @(x)isnumeric(x)&&x>0);
parse(p,varargin{:});
inputSpatial = p.Results.inputSpatial;
outputSpatial = p.Results.outputSpatial;
fullResX = p.Results.fullResX;
hexaSampling = p.Results.hexaSampling;
disableGPU = p.Results.disableGPU;
nbImgGPU = p.Results.nbImgGPU;
nbRefocusGPU = p.Results.nbRefocusGPU;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~exist('padParams','var')), padParams = BorderParams();end
useGPU = parallel.gpu.GPUDevice.isAvailable && ~disableGPU;


nChan = size(Views,3);
numViews = size(Views,4);

if(inputSpatial || isempty(fullResX) )
    fullResX = size(Views,2);
end
if(inputSpatial)%only perform padding if input is in spatial domain.
    fullResX = fullResX + padParams.L + padParams.R;
    fullResY = size(Views,1) + padParams.T + padParams.B ;
else
    fullResY = size(Views,1);
end
[wx,wy,xC,~] = GenerateFrequencyGrid(fullResX,fullResY);

if(useGPU)
    U=gpuArray(reshape(U,1,1,1,1,numViews));
    V=gpuArray(reshape(V,1,1,1,1,numViews));
    dRange = gpuArray(reshape(dRange,1,1,1,[]));
    wx = gpuArray(wx(:,1:xC));
    wy = gpuArray(wy(:,1:xC));
else
    wx = wx(:,1:xC);
    wy = wy(:,1:xC);
end

%% Transform input views (include padding and windowing).
if(inputSpatial)
    Views = fftImgs(Views, padParams, 'hexaSampling', hexaSampling, 'disableGPU', disableGPU, 'nbImgGPU', nbImgGPU);
end
Views = Views(:,1:xC,:,:);

%% Perform shift sum refocusing in Fourier domain
stackSize = length(dRange);
I = zeros([fullResY xC nChan stackSize],'single');

if(useGPU)
    Views=reshape(Views,fullResY,xC,nChan,1,[]);
    for batchStart=1:nbRefocusGPU:stackSize
        id_Disps = batchStart:min(batchStart+nbRefocusGPU-1,stackSize);
        for l=1:nbImgGPU:numViews
            idViews = l:min(l+nbImgGPU-1,numViews);
            Imgs  = gpuArray(Views(:,:,:,:,idViews));
            I(:,:,:,id_Disps) = I(:,:,:,id_Disps) + gather(sum( bsxfun(@times, Imgs, arrayfun(@renderShiftWeights,wx,wy,U(:,:,:,:,idViews),V(:,:,:,:,idViews),dRange(id_Disps))),5));
        end
    end
else
    for id_d=1:stackSize
        d = dRange(id_d);
        for l=1:numViews
            Shift = exp(-2i*pi*d*(wx*U(l) + wy*V(l)));
            for ch=1:nChan
                I(:,:,ch,id_d) = I(:,:,ch,id_d) + Views(:,:,ch,l).*Shift;
            end
        end
    end
end
I = I./numViews;

%% Inverse transform the result (include cropping padded borders if output is in spatial domain).
if(outputSpatial)
    I = ifftImgs(I,padParams,'fullResX', fullResX, 'disableGPU', disableGPU, 'nbImgGPU', nbImgGPU);
end

end


function out = renderShiftWeights(wx,wy, u, v,disp)
    out = exp(-2i*pi*disp*(wx*u+wy*v));
end
