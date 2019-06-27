% Gradient-descent based Calibration from a relaxed version of the FDL
% model with no Rank-1 constraint on the parameter matrix (containing the
% shift of each layer to reconstruct each view).
% Note: FDL constructed from this relaxed model are only suitable to
% reconstruct images at the same view positions as the input views (e.g.
% for denoising).
%
%Inputs:
% - ViewsFFT : Fourier Transform of input views with dimensions : (#views, #color channels, Y resolution, X resolution).
% - wx,wy : Horizontal and vertical frequency values associated to each frequency component in ViewsFFT. Dimensions = (Y resolution, X resolution).
% - numDisp: Number of disparity values (i.e. number of layers).
% - lambda : Regularization parameter (typical value = 1).
% - PxInit : (Optional) Initialization of the matrix of parameters for the horizontal dimension.
% - PyInit : (Optional) Initialization of the matrix of parameters for the vertical dimension.
%
%Outputs:
% - Px : Matrix of parameters for the horizontal dimension (without rank 1 constraint).
% - Py : Matrix of parameters for the vertical dimension (without rank 1 constraint).
% - U, V, D : Lists of view positions and disparity values derived from rank 1 approximation of [Px;Py].
% - sigPx,sigPy : Standard deviation of each parameter in Px,Py estimated by observing the parameters'
% variation over a given number of iterations (numIterSigError) after the convergence criterion is reached.
% - residN2 : Squared residual at each iteration.
% - gradN2 : Squared gradient norm for D, U and V at each iteration.

function [Px,Py,U,V,D,sigPx,sigPy,residN2,gradN2,gX,gY] = CalibrateFDL_NoR1_cpu(ViewsFFT,wx,wy,numDisp,lambda,PxInit,PyInit)

ViewsFFT = single(ViewsFFT);

nChan = size(ViewsFFT,2);
if(nChan==3)
    ViewsFFT = 0.2126*ViewsFFT(:,1,:,:) + 0.7152*ViewsFFT(:,2,:,:) + .0722*ViewsFFT(:,3,:,:);
    nChan=1;
end

xC = ceil((size(ViewsFFT,4)+1)/2);
yC = ceil((size(ViewsFFT,3)+1)/2);

numInitViews = size(ViewsFFT,1);
numFreqX = size(ViewsFFT,4);
numFreqY = size(ViewsFFT,3);
numFreqProc = numFreqY*(xC-1)+yC -1;%first half of the spectrum, omitting the frequency (0,0).

wx = single(wx(:,1:xC));
wx = permute(wx(:),[3 2 1]);
wy = single(wy(:,1:xC));
wy = permute(wy(:),[3 2 1]);

%Adjustment of the regularization parameter to keep the results consistent with
%repect to number of input views and number of layers.
lambda = lambda*numDisp*numInitViews/4000;

%Matrix for 2nd order layer regularization
L2=full(gallery('tridiag',numDisp,1,-2,1));
L2(2:end+1,:)=L2;L2(1,1:2)=[1 0];L2(end+1,end-1:end)=[0 1];
Reg = single(lambda*(L2'*L2));%2nd order regularization term


%Initialization
if(exist('PxInit','var') && size(PxInit,1)==numInitViews && size(PxInit,2)==numDisp)
    Px = single(PxInit);
else
    Px = zeros(numInitViews,numDisp,'single');
end
if(exist('PyInit','var') && size(PyInit,1)==numInitViews && size(PyInit,2)==numDisp)
    Py = single(PyInit);
else
    Py = zeros(numInitViews,numDisp,'single');
end

gradPx = zeros(numInitViews,numDisp);
gradPy = zeros(numInitViews,numDisp);

rng('default');rng(0);

%Number of frequencies used in a mini-batch.
freqBlockSz = 4096;%16384;%

%Variable used to determine the order of magnitude of a "small" value with respect to the gradients (used for "attenuated" normalization).
EpsgradP  = 1e-14 * 4*pi*2*numDisp * freqBlockSz * sum(abs(ViewsFFT(:,:,yC,xC)).^2);

%Parameters for stopping criterion:
EvalPeriod = 1;    %Evaluate residual (and stop criterion) every 'EvalPeriod' iterations.
numIterSmooth = 10;%Number of iterations for residual smoothing over iterations (higher => more reliable stopping criterion).
stopTh = 1e-4;     %Threshold on the relative change of the smoothed residual, for the stopping criterion.
freqEvalIds = randperm(numFreqProc,freqBlockSz); %Pre-selected frequencies to evaluate the residual for stopping criterion.
rng('default');rng(0);

%Parameters for estimation of the mean and standard deviation of the parameters.
numIterSigError=10;%Number of iteration after convergence to estimate average and standard deviation of the estimated parameters.
EvaluatingSigError = false;
sumPx=zeros(size(Px)); sigPx=zeros(size(Px));
sumPy=zeros(size(Py)); sigPy=zeros(size(Py));

%Main Loop
freqBlockSz=min(freqBlockSz,numFreqProc);
numIterMax=500;

iter=1;
resid = zeros(numInitViews, nChan, freqBlockSz);
Z = zeros(numDisp, nChan, freqBlockSz);
while(iter <= numIterMax)
    freqIds = randperm(numFreqProc,freqBlockSz);
    
%%%%%%%%%%%%%%%%%%%%%%%
%  Compute Gradients  %
%%%%%%%%%%%%%%%%%%%%%%%

    A = exp( 2i*pi*(bsxfun(@times,Px,wx(freqIds)) + bsxfun(@times,Py,wy(freqIds))) );
    for i = 1:freqBlockSz
        freqId = freqIds(i);
        Z(:,:,i) = (A(:,:,i)'*A(:,:,i) + Reg) \ ( A(:,:,i)' * ViewsFFT(:,:,freqId));
        resid(:,:,i) = A(:,:,i) * Z(:,:,i) - ViewsFFT(:,:,freqId);
    end
    
    for k=1:numDisp
        gradTmp = 4*pi*imag(bsxfun(@times,...
                            conj(bsxfun(@times,Z(k,:,:),A(:,k,:))),...
                            resid));
        gradPx(:,k) = sum(sum(bsxfun(@times, wx(freqIds), gradTmp), 2), 3);
        gradPy(:,k) = sum(sum(bsxfun(@times, wy(freqIds), gradTmp), 2), 3);
    end

    gradN2(iter,1) = sum(gradPx(:).^2);
    gradN2(iter,2) = sum(gradPy(:).^2);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Evaluate Residual and stopping criterion          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%residN2(iter) = sum(abs(resid(:)).^2);
%%{
if(mod(iter,EvalPeriod)==0 && ~EvaluatingSigError)
    idEval=iter/EvalPeriod;
    A = exp( 2i*pi*(bsxfun(@times,Px,wx(freqEvalIds)) + bsxfun(@times,Py,wy(freqEvalIds))) );
    for i=1:freqBlockSz
        freqId = freqEvalIds(i);
        Z(:,:,i) = (A(:,:,i)'*A(:,:,i) + Reg) \ ( A(:,:,i)' * ViewsFFT(:,:,freqId));
        resid(:,:,i) = A(:,:,i) * Z(:,:,i) - ViewsFFT(:,:,freqId);
    end
    
    residN2(idEval,1) = sum(abs(resid(:)).^2);
    if(idEval>2*numIterSmooth-1),residN2(idEval,2) = mean(residN2(idEval-numIterSmooth+1:idEval,1));end
    if(idEval>2*numIterSmooth-1 && abs(mean(residN2(idEval-numIterSmooth+1:idEval,1))-mean(residN2(idEval-2*numIterSmooth+1:idEval-numIterSmooth,1))) < mean(residN2(idEval-numIterSmooth+1:idEval,1))*stopTh)
        disp(['Criterion met after ' num2str(iter) ' iterations.']);
        numIterMax = iter + numIterSigError;
        EvaluatingSigError = true;
        disp(['Estimating mean and standard deviation of the parameters (' num2str(numIterSigError) ' more iterations) ...']);
    end
end
if(iter == numIterMax)
    if(EvaluatingSigError)
        Px = sumPx/numIterSigError; Py = sumPy/numIterSigError;
        sigPx = sqrt(max(0,(sigPx-sumPx.*sumPx/numIterSigError))/numIterSigError);
        sigPy = sqrt(max(0,(sigPy-sumPy.*sumPy/numIterSigError))/numIterSigError);
        break;
    else
        disp(['No convergence after the maximum number of iterations (' num2str(numIterMax) ').']);
        numIterMax = iter + numIterSigError;
        EvaluatingSigError = true;
        disp(['Estimating mean and standard deviation of the parameters (' num2str(numIterSigError) ' more iterations) ...']);
    end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%
%     Line search     %
%%%%%%%%%%%%%%%%%%%%%%%
    stepPx = 1/(EpsgradP+sqrt(sum(gradN2(iter,:))));
    stepPy = stepPx;

%%%%%%%%%%%%%%%%%%%%%%%
%       Update        %
%%%%%%%%%%%%%%%%%%%%%%%

    Px = Px - stepPx*gradPx;
    Py = Py - stepPy*gradPy;

    [~,~,Svr] = svd([Px;Py],'econ');
    [~, Ids] = sort(Svr(:,1));
    Px = Px(:,Ids);
    Py = Py(:,Ids);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute standard deviation and average parameters with additional iterations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(EvaluatingSigError)
        sumPx = sumPx + Px; sigPx = sigPx + Px.*Px;
        sumPy = sumPy + Py; sigPy = sigPy + Py.*Py;
    end


iter=iter+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Px = gather(Px);
Py = gather(Py);

gradN2 = gather(gradN2);
residN2 = gather(residN2);

gX = gather(stepPx*gradPx);
gY = gather(stepPy*gradPy);

% Rank 1 approximation of [Px;Py] to find positions U,V and
% disparities D (with rank 1 approx: [Px;Py] = [U;V]*D).
[Svl,S,Svr] = svd([Px;Py],'econ');
U = Svl(1:size(Px,1),1)*S(1,1);
V = Svl(1+size(Px,1):end,1)*S(1,1);
D = (Svr(:,1))';

end