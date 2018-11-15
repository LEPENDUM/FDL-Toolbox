LFDir = [fileparts(mfilename('fullpath')) '/'];
gpuInit = true;%true -> use gpu for initilisation (Windowing + Fourier Transform of input views).

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
numImages = nU*nV;


InitViewIds=[1:numImages];%Selected views for the input of the FDL construction.
numInitViews = length(InitViewIds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parameters for FDL calibration + construction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambdaCalib = 1;
lambdaConstruct = .1;
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


%% Prepare variables
fullResX = padSizeX*2+imgSize(2);
fullResY = padSizeY*2+imgSize(1);

%frequency grid:
[wx,wy,xC,yC] = GenerateFrequencyGrid(fullResX,fullResY);

%window profile:
if(Windowing)
    Window = single(GenerateWindowIntensity(fullResX,fullResY,padSizeX,padSizeY, 'hann'));
    if(gpuInit), Window = gpuArray(Window);end
end



%% Pre-process input data (Padding/Windowing + Fourier Tramsform )
ViewsFFT = zeros([numInitViews,nChan,fullResY,fullResX],'single');
tic;fprintf('Computing Fourier Transform of input views...');
for l=1:numInitViews
    idView = InitViewIds(l);
    
    if(gpuInit)
        orgPad = gpuArray(padarray(single(LF(:,:,:,idView)),[padSizeY padSizeX],padType));
    else
        orgPad = padarray(single(LF(:,:,:,idView)),[padSizeY padSizeX],padType);
    end
    if(Windowing)
        orgPad = bsxfun(@times, orgPad, Window);
    end
    
    if(gpuInit)
        ViewsFFT(l,:,:,:) = gather(permute(fftshift(fftshift(fft2(orgPad),1),2),[3 1 2]));
    else
        for ch=1:nChan, ViewsFFT(l,ch,:,:) = fftshift(fft2(orgPad(:,:,ch))); end
    end
end
t=toc;fprintf(['\b\b\b (done in ' num2str(t) 's)\n']);

clear Window orgPad


%% FDL Calibration
tic;fprintf('FDL Calibration...');
[Px,Py,U,V,D,~,~,~,~,~,residN2,gradN2]=CalibrateFDL_UVD_gpu(ViewsFFT, wx, wy, numLayers, lambdaCalib);
[U,V,D]=ScaleParams(U,V,D,nU,nV); %Scaling and centering of the parameters to match approximately an integer grid of views (this step has no effect on the final images).
t=toc;fprintf(['\b(' num2str(t) 's)\n']);

%% FDL Construction
tic;fprintf('FDL Construction...');
FDL = ComputeFDL_gpu(ViewsFFT, wx, wy, U, V, D, lambdaConstruct);
t=toc;fprintf(['(' num2str(t) 's)\n']);

%% Start Rendering Application
RenderAppMain(FDL, [fullResY,fullResX], [padSizeX, padSizeX, padSizeY, padSizeY], D,U,V, [],{linearizeInput,gammaOffset});
