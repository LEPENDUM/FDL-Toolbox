% Gradient-descent based Calibration from the FDL model.
% The function determines the positions of the input views on the camera
% planes and a set of disparity values needed for FDL construction.
%
%Inputs:
% - ViewsFFT : Fourier Transform of input views with dimensions: 1.Views, 2.Color channels, 3.Vertical axis (Y), 4.Horizontal axis (X).
%
% - wx : Horizontal frequency values associated to each frequency component in ViewsFFT. size = (Y resolution, X resolution).
%
% - wy : Vertical frequency values associated to each frequency component in ViewsFFT. size = (Y resolution, X resolution).
%
% - numDisp: Number of disparity values (i.e. number of layers).
%
% - lambda : Regularization parameter (typical value = 1).
%
% - uPosInit : (Optional) Initialisation for horizontal positions of the views (length(uPosInit)= #views).
%              The number of elements in the 2 first dimensions combined must be the number of views.
%              If the option SeparateChannels is active, the number of elements in the 3rd dimension may be either 1 or the number of channels.
%
% - vPosInit : (Optional) Initialisation for vertical positions of the views (length(vPosInit)= #views).
%              The number of elements in the 2 first dimensions combined must be the number of views.
%              If the option SeparateChannels is active, the number of elements in the 3rd dimension may be either 1 or the number of channels.
%
% - DispsInit : (Optional) Initialisation for the layer's disparities (length(DispsInit)=numDisp).
%              The number of elements in the 2 first dimensions combined must be the number of layers.
%              If the option SeparateChannels is active, the number of elements in the 3rd dimension may be either 1 or the number of channels.
%
% - varargin: Additional optional arguments given as (Name,Value) pairs:
%   - 'SeparateChannels' (default=false) : 
%       true->Compute separate calibration parameters for each color channel.
%       false->perform calibration on luma (if RGB components are given).
%   - 'fixed_D' (default=false) :
%       Set to true to prevent update of the disparity values (requires the disparity values to be given as input in the DispsInit argument).
%   - 'fixed_UV' (default=false) :
%       Set to true to prevent update of the view coordinates (requires the view coordinates to be given as input in the uPosInit and vPosInit arguments).
%   - 'UseSignFlip' (default=true) :
%       Flip the sign of view coordinates in one dimension at early iterations to reduce the risk of falling in local minimum (only applies if fixed_UV is false).
%   - 'numIterSignEval' (default=15) :
%       Number of initial iterations for sign flip evaluation (the sign giving lowest error is decided afterwards and kept for the remaining iterations).
%   - 'numIterLowFreqs' (default=20) :
%       Number of initial iterations for which only low frequencies are selected for gradient computation and parameter update.
%   - 'lowFreqThreshold' (default=0.1) :
%       Threshold for the low frequency selection between 0 (no frequency) and 0.5 (all the frequencies).
%   - 'mainChan' (default=0) :
%       Index of the main channel (for the stopping criterion evaluation, only the calibration parameters from the main channel are used).
%       Set to 0 for the automatic default choice: 1 if only one component or if SeparateChannels is false / 2 otherwise (i.e. green component of RGB input).
%   - 'stopTh' (default=1e-3) : Threshold on the relative change of the smoothed residual, for the stopping criterion.
%   - 'numIterMax' (default=500) : Maximum number of iterations.
%   - 'numIterSigError' (default=10) :
%       Number of iterations after the stopping criterion is met to estimate standard deviation of the estimated parameters.
%       The final parameters returned are averaged over these last iterations.
%   - 'numIterSmooth' (default=10) :
%       Number of iterations for residual smoothing over iterations (for more robust stopping criterion).
%   - 'EvalPeriod' (default=1) :
%       Evaluate residual (and stopping criterion) every 'EvalPeriod' iterations.
%   - 'freqBlockSz' (default=4096)
%       Number of frequencies used in a mini-batch (only if it is lower or equal to the number of frequencies avaialble in the data).
%
%
%Outputs:
% - Px : Rank 1 matrix of parameters for the horizontal dimension.
%
% - Py : Rank 1 matrix of parameters for the vertical dimension.
%
% - U, V, D : Lists of view positions and disparity values.
% The dimensions of U and V are (#views, 1, #channels) if the option SeparateChannels is active, and (#views, 1) otherwise.
% The dimensions of D are (1, #layers, #channels) if the option SeparateChannels is active, and (1, #layers) otherwise.
%
% - sigPx,sigPy,sigU,sigV,sigD : Standard deviation of each parameter in Px,Py,U,V,D estimated by observing the 
% parameters' variation over a given number of iterations (numIterSigError) after the convergence criterion is reached.
%
% - residN2 : Squared residual at each iteration.
%
% - gradN2 : Squared gradient norm for D, U and V at each iteration.

function [Px,Py,U,V,D,sigPx,sigPy,sigU,sigV,sigD,residN2,gradN2] = CalibrateFDL_UVD_cpu(ViewsFFT,wx,wy,numDisp,lambda, uPosInit,vPosInit,DispsInit,varargin)

p = inputParser;
addParameter (p,'SeparateChannels',     false,  @islogical);
addParameter (p,'fixed_D',              false,  @islogical);
addParameter (p,'fixed_UV',             false,  @islogical);
addParameter (p,'UseSignFlip',          true,   @islogical);
addParameter (p,'numIterSignEval',      15,     @isnumeric);
addParameter (p,'numIterLowFreqs',      20,     @isnumeric);
addParameter (p,'lowFreqThreshold',     0.1,    @isnumeric);
addParameter (p,'numIterSigError',       10,    @isnumeric);

addParameter (p,'mainChan',             0,      @isnumeric);
addParameter (p,'stopTh',               1e-3,   @isnumeric);
addParameter (p,'numIterMax',           500,    @isnumeric);
addParameter (p,'numIterSmooth',        10,     @isnumeric);
addParameter (p,'EvalPeriod',           1,      @isnumeric);
addParameter (p,'freqBlockSz',          4096,   @isnumeric);


parse(p,varargin{:});

numIterMax = p.Results.numIterMax;
EvalPeriod = p.Results.EvalPeriod;
numIterSmooth = p.Results.numIterSmooth;
stopTh = p.Results.stopTh;
mainChan = p.Results.mainChan;
UseSignFlip = p.Results.UseSignFlip;
numIterSignEval = p.Results.numIterSignEval;
numIterLowFreqs = p.Results.numIterLowFreqs;
lowFreqThreshold = p.Results.lowFreqThreshold;
freqBlockSz = p.Results.freqBlockSz;
fixed_D = p.Results.fixed_D;
fixed_UV = p.Results.fixed_UV;
SeparateChannels = p.Results.SeparateChannels;
numIterSigError = p.Results.numIterSigError;

if(fixed_UV), UseSignFlip=false;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ViewsFFT = single(ViewsFFT);

nChan = size(ViewsFFT,2);
if(~SeparateChannels)
    if(nChan==3)
        ViewsFFT = 0.2126*ViewsFFT(:,1,:,:) + 0.7152*ViewsFFT(:,2,:,:) + .0722*ViewsFFT(:,3,:,:);
        nChan=1;
    end
    if(nChan~=1), error('When the ''SeparateChannels'' option is not active, the number of channels must be either 1 or 3 (RGB components).');end
end
if(mainChan==0)
    mainChan=min(2,nChan);
elseif(mainChan>nChan)
    if(nChan>1),error('The selected mainChan parameter should be lower or equal to the number of channels in the input data.');
    else, mainChan=1;
    end
end

xC = ceil((size(ViewsFFT,4)+1)/2);
yC = ceil((size(ViewsFFT,3)+1)/2);

numInitViews = size(ViewsFFT,1);
numFreqX = size(ViewsFFT,4);
numFreqY = size(ViewsFFT,3);
numFreqProc = numFreqY*(xC-1)+yC -1;%First half of the spectrum, omitting the frequency (0,0).

wx = single(wx(:,1:xC));
wx = permute(wx(:),[4 3 2 1]);
wy = single(wy(:,1:xC));
wy = permute(wy(:),[4 3 2 1]);

%Adjustment of the regularization parameter to keep the results consistent with
%repect to number of input views and number of layers.
lambda = lambda*numDisp*numInitViews/4;

%Matrix for 2nd order layer regularization
if(numDisp>2)
    L2=full(gallery('tridiag',numDisp,1,-2,1));
    L2(2:end+1,:)=L2;L2(1,1:2)=[1 0];L2(end+1,end-1:end)=[0 1];
elseif(numDisp==2)
    L2 = [-2 1; 1 -2];
elseif(numDisp==1)
    L2 = 1;
else
    error('Parameter ''numDisp'' should be a strictly positive integer.');
end
Reg = single(lambda*(L2'*L2));


%Initializations
if(exist('uPosInit','var') && ~isempty(uPosInit))
    if( size(uPosInit,1)*size(uPosInit,2)==numInitViews && (size(uPosInit,3)==nChan || size(uPosInit,3)==1) )
        U = bsxfun(@times, reshape(single(uPosInit),numInitViews,1,size(uPosInit,3)), ones(1,1,nChan) );
    else
        error('Incorrect dimensions of input parameter uPosInit.');
    end
elseif(fixed_UV)
    error('The uPosInit and vPosInit argument are required when the ''fixed_UV'' option is active.');
else
    U = zeros(numInitViews,1,nChan,'single');
end
if(exist('vPosInit','var') && ~isempty(uPosInit))
    if( size(vPosInit,1)*size(vPosInit,2)==numInitViews && (size(vPosInit,3)==nChan || size(vPosInit,3)==1) )
        V = bsxfun(@times, reshape(single(vPosInit),numInitViews,1,size(vPosInit,3)), ones(1,1,nChan) );
    else
        error('Incorrect dimensions of input parameter vPosInit.');
    end
elseif(fixed_UV)
    error('The uPosInit and vPosInit argument are required when the ''fixed_UV'' option is active.');
else
    V = zeros(numInitViews,1,nChan,'single');
end
if(exist('DispsInit','var') && ~isempty(DispsInit))
    if( size(DispsInit,1)*size(DispsInit,2)==numDisp && (size(DispsInit,3)==nChan || size(DispsInit,3)==1) )
        D = bsxfun(@times, reshape(single(DispsInit),1,numDisp,size(DispsInit,3)), ones(1,1,nChan) );
    else
        error('Incorrect dimensions of input parameter DispsInit.');
    end
elseif(fixed_D)
    error('The DispsInit argument is required when the ''fixed_D'' option is active.');    
else
    D = repmat(linspace(0,10,numDisp),1,1,nChan);
end
if(~fixed_D),  gradD = zeros(size(D),'single');end
if(~fixed_UV), gradU = zeros(size(U),'single');
               gradV = zeros(size(V),'single');
end
gradN2=single([]);


%Adjust mini-batch size if the number of frequencies in the input data is too small.
freqBlockSz = min(freqBlockSz,numFreqProc);

%Variables used to determine the order of magnitude of a "small" value with respect to the gradients (used for "attenuated normalization").
EpsgradD  = 1e-13 * 4*pi*2*numDisp * freqBlockSz * sum(sum(abs(ViewsFFT(:,:,yC,xC)).^2));
EpsgradUV = EpsgradD;

%Pre-selected frequency sets:
    %Low frequencies for update at first iterations.
freqsProcInit = find((wx(1:numFreqProc).^2+wy(1:numFreqProc).^2) < lowFreqThreshold);
numFreqProcInit = numel(freqsProcInit);
    %Frequencies to evaluate the residual for stopping criterion.
rng('default');rng(0);
freqEvalIds = randperm(numFreqProc,freqBlockSz);

%Mean and standard deviation of the parameters at iterations after reaching the stopping criterion.
sumPx=zeros(length(U),length(D),nChan); sigPx=zeros(length(U),length(D),nChan);
sumPy=zeros(length(U),length(D),nChan); sigPy=zeros(length(U),length(D),nChan);
sumU=zeros(size(U));   sigU=zeros(size(U));
sumV=zeros(size(V));   sigV=zeros(size(V));
sumD=zeros(size(D));   sigD=zeros(size(D));

%Main Loop
iter=1;
finishedMainIter=false; EvaluatingSigError = false;
resid = zeros(numInitViews, 1, nChan, freqBlockSz, 'single');
Z = zeros(numDisp, 1, nChan, freqBlockSz, 'single');
rng('default');rng(0);
while(iter <= numIterMax)
    if(iter<=numIterLowFreqs)
        freqIds = freqsProcInit(randperm(numFreqProcInit,min(freqBlockSz,numFreqProcInit)));
    else
        freqIds = randperm(numFreqProc,freqBlockSz);
    end
    
    %Test sign flipping at first iterations
    if(UseSignFlip)
        if(iter>1 && iter<=numIterSignEval)
            V=-V;
        end
        %Decide the final sign from residuals at previous iterations.
        if(iter==numIterSignEval+1)
            residN2Smooth = (residN2(2:numIterSignEval-2,1)+residN2(3:numIterSignEval-1,1)+residN2(4:numIterSignEval,1))/3; %smooth residual from 3 latest iterations
            residN2Smooth = residN2(3:numIterSignEval-1,1) - residN2Smooth;                                                 %Difference between true and smoothed residual
            if(mean(residN2Smooth(end:-2:1)) < mean(residN2Smooth(end-1:-2:1)))
                V=-V;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%
%  Compute Gradients  %
%%%%%%%%%%%%%%%%%%%%%%%

    A = exp( 2i*pi*(bsxfun(@times,bsxfun(@times,U,D),wx(freqIds)) + bsxfun(@times,bsxfun(@times,V,D),wy(freqIds))) );
    
    for i = 1:freqBlockSz
        freqId = freqIds(i);
        for ch=1:nChan
            Z(:,:,ch,i) = (A(:,:,ch,i)'*A(:,:,ch,i) + Reg) \ ( A(:,:,ch,i)' * ViewsFFT(:,ch,freqId));
            resid(:,:,ch,i) = A(:,:,ch,i) * Z(:,:,ch,i) - ViewsFFT(:,ch,freqId);
        end
    end
    
    gradTmp = 4*pi*imag(bsxfun(@times,...
            conj(bsxfun(@times,permute(Z,[2 1 3 4]),A)),...
            resid));
    clear A
    
    gradPx = sum(bsxfun(@times, wx(freqIds), gradTmp), 4);
    gradPy = sum(bsxfun(@times, wy(freqIds), gradTmp), 4);
    clear gradTmp
    
    if(~fixed_D)
        gradD = sum(bsxfun(@times,gradPx,U) + bsxfun(@times,gradPy,V),1) +.3*gradD;
        gradN2(iter,1) = sum(gradD(:).^2);
    end
    if(~fixed_UV)
        gradU = sum(bsxfun(@times,gradPx,D),2) +.3*gradU;
        gradV = sum(bsxfun(@times,gradPy,D),2) +.3*gradV;
        gradN2(iter,2) = sum(gradU(:).^2);
        gradN2(iter,3) = sum(gradV(:).^2);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Evaluate Residual and stopping criterion          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mod(iter,EvalPeriod)==0 && ~EvaluatingSigError)
    idEval = iter/EvalPeriod;
    A = exp( 2i*pi*(bsxfun(@times,U(:,:,mainChan)*D(:,:,mainChan),wx(freqEvalIds)) + bsxfun(@times,V(:,:,mainChan)*D(:,:,mainChan),wy(freqEvalIds))) );
    for i = 1:freqBlockSz
        freqId = freqEvalIds(i);
        Z(:,:,:,i) = (A(:,:,1,i)'*A(:,:,1,i) + Reg) \ ( A(:,:,1,i)' * ViewsFFT(:,:,freqId) );
        resid(:,:,:,i) = A(:,:,1,i) * reshape(Z(:,:,:,i),[],nChan) - ViewsFFT(:,:,freqId);
    end
    residN2(idEval,1) = sum(abs(resid(:)).^2);
    if(idEval>2*numIterSmooth-1),residN2(idEval,2) = mean(residN2(idEval-numIterSmooth+1:idEval,1));end
    if(idEval>2*numIterSmooth-1 && abs(mean(residN2(idEval-numIterSmooth+1:idEval,1))-mean(residN2(idEval-2*numIterSmooth+1:idEval-numIterSmooth,1))) < mean(residN2(idEval-numIterSmooth+1:idEval,1))*stopTh)
        disp(['Criterion met after ' num2str(iter) ' iterations.']);
        numIterMax = iter + numIterSigError;
        finishedMainIter = true;
        EvaluatingSigError = numIterSigError > 0;
        if(EvaluatingSigError)
            disp(['Estimating mean and standard deviation of the parameters (' num2str(numIterSigError) ' more iterations) ...']);
        end
    end
end
if(iter == numIterMax)
    if(EvaluatingSigError)
        if(~fixed_UV),  U = sumU/numIterSigError; V = sumV/numIterSigError;
                        sigU = sqrt(max(0,(sigU-sumU.*sumU/numIterSigError))/numIterSigError);
                        sigV = sqrt(max(0,(sigV-sumV.*sumV/numIterSigError))/numIterSigError);
        end
        if(~fixed_D),   D = sumD/numIterSigError;
                        sigD = sqrt(max(0,(sigD-sumD.*sumD/numIterSigError))/numIterSigError);
        end
        Px = sumPx/numIterSigError; Py = sumPy/numIterSigError;
        sigPx = sqrt(max(0,(sigPx-sumPx.*sumPx/numIterSigError))/numIterSigError);
        sigPy = sqrt(max(0,(sigPy-sumPy.*sumPy/numIterSigError))/numIterSigError);
        break;
    elseif(~finishedMainIter)
        disp(['No convergence after the maximum number of iterations (' num2str(numIterMax) ').']);
        numIterMax = iter + numIterSigError;
        finishedMainIter = true;
        EvaluatingSigError = numIterSigError > 0;
        if(EvaluatingSigError)
            disp(['Estimating mean and standard deviation of the parameters (' num2str(numIterSigError) ' more iterations) ...']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Line search and update     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(~fixed_D)
        %Line search
        stepD = .2/( EpsgradD + sqrt(gradN2(iter,1)));
        %Update
        D = D - stepD*gradD;
    end
    if(~fixed_UV)
        %Line search
        stepU = .2/( EpsgradUV + sqrt(sum(gradN2(iter,2:3))));
        stepV = stepU;
        %Update
        U = U - stepU*gradU;
        V = V - stepV*gradV;
        if(~fixed_D)
        %Normalize U,V and D => the normalized solution is strictly equivlalent w.r.t.
        %the objective function, but the normalization improves the convergence.
            normUV = max(max(U(:))-min(U(:)),max(V(:))-min(V(:)))/2;
            D = sort(D) * normUV;
            U = U / normUV;
            V = V / normUV;
        end
    elseif(~fixed_D)
        D = sort(D);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute standard deviation and average parameters with additional iterations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(EvaluatingSigError)
        Px = bsxfun(@times,U,D); Py = bsxfun(@times,V,D);
        sumPx = sumPx + Px; sigPx = sigPx + Px.*Px;
        sumPy = sumPy + Py; sigPy = sigPy + Py.*Py; 
        if(~fixed_UV), sumU = sumU + U;    sigU  = sigU + U.*U;
                       sumV = sumV + V;    sigV  = sigV + V.*V;
        end
        if(~fixed_D),  sumD = sumD + D;    sigD  = sigD + D.*D;end
    end


iter=iter+1;
end

if(~EvaluatingSigError)
    Px = bsxfun(@times,U,D);
    Py = bsxfun(@times,V,D);
end

disp('done!');
end
