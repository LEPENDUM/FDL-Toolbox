% Gradient-descent based Calibration from the FDL model.
% The function determines the positions of the input views on the camera
% planes and a set of disparity values needed for FDL construction.
%
%Inputs:
% - ViewsFFT : Fourier Transform of input views with dimensions : (#views, #color channels, Y resolution, X resolution).
% - wx,wy : Horizontal and vertical frequency values associated to each frequency component in ViewsFFT. Dimensions = (Y resolution, X resolution).
% - numDisp: Number of disparity values (i.e. number of layers).
% - lambda : Regularization parameter (typical value = 1).
% - uPosInit : (Optional) Initialisation for horizontal positions of the views (length(uPosInit)= #views).
% - vPosInit : (Optional) Initialisation for vertical positions of the views (length(vPosInit)= #views).
% - DispsInit : (Optional) Initialisation for the layer's disparities (length(DispsInit)=numDisp).
%
%Outputs:
% - Px : Rank 1 matrix of parameters for the horizontal dimension.
% - Py : Rank 1 matrix of parameters for the vertical dimension.
% - U, V, D : Lists of view positions and disparity values.
% - sigPx,sigPy,sigU,sigV,sigD : Standard deviation of each parameter in Px,Py,U,V,D estimated by observing the 
% parameters' variation over a given number of iterations (numIterSigError) after the convergence criterion is reached.
% - residN2 : Squared residual at each iteration.
% - gradN2 : Squared gradient norm for D, U and V at each iteration.

function [Px,Py,U,V,D,sigPx,sigPy,sigU,sigV,sigD,residN2,gradN2] = CalibrateFDL_UVD_cpu(ViewsFFT,wx,wy,numDisp,lambda, uPosInit,vPosInit,DispsInit)

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
numFreqProc = numFreqY*(xC-1)+yC -1;%First half of the spectrum, omitting the frequency (0,0).

wx = single(wx(:,1:xC));
wx = permute(wx(:),[3 2 1]);
wy = single(wy(:,1:xC));
wy = permute(wy(:),[3 2 1]);

%Adjustment of the regularization parameter to keep the results consistent with
%repect to number of input views and number of layers.
lambda = lambda*numDisp*numInitViews/4;

%Matrix for 2nd order layer regularization
L2=full(gallery('tridiag',numDisp,1,-2,1));
L2(2:end+1,:)=L2;L2(1,1:2)=[1 0];L2(end+1,end-1:end)=[0 1];
Reg = single(lambda*(L2'*L2));


%Initialization
gradPx = zeros(numInitViews,numDisp);
gradPy = zeros(numInitViews,numDisp);

rng('default');rng(0);

if(exist('uPosInit','var') && numel(uPosInit)==numInitViews)
    U = single(uPosInit(:));
else
    U = zeros(numInitViews,1,'single');
end
if(exist('vPosInit','var') && numel(vPosInit)==numInitViews)
    V = single(vPosInit(:));
else
    V = zeros(numInitViews,1,'single');
end
if(exist('DispsInit','var') && numel(DispsInit)==numDisp)
    D = single(DispsInit(:)');
else
    D = linspace(0,10,numDisp);
end
gradD = zeros(size(D),'single');
gradU = zeros(size(U),'single');
gradV = zeros(size(V),'single');




%Number of frequencies used in a mini-batch.
freqBlockSz = 4096;

%Variables used to determine the order of magnitude of a "small" value with respect to the gradients (used for "attenuated" normalization).
EpsgradD  = 1e-13 * 4*pi*2*numDisp * freqBlockSz * sum(abs(ViewsFFT(:,:,yC,xC)).^2);
EpsgradUV = EpsgradD;

%Parameters for stopping criterion:
EvalPeriod = 1;    %Evaluate residual (and stop criterion) every 'EvalPeriod' iterations.
numIterSmooth = 10;%Number of iterations for residual smoothing over iterations (higher => more reliable stopping criterion).
stopTh = 1e-3;     %Threshold on the relative change of the smoothed residual, for the stopping criterion.
freqEvalIds = randperm(numFreqProc,freqBlockSz); %Pre-selected frequencies to evaluate the residual for stopping criterion.
rng('default');rng(0);

%Parameters for estimation of the mean and standard deviation of the parameters.
numIterSigError=10;%Number of iteration after convergence to estimate average and standard deviation of the estimated parameters.
EvaluatingSigError = false;
sumPx=zeros(length(U),length(D)); sigPx=zeros(length(U),length(D));
sumPy=zeros(length(U),length(D)); sigPy=zeros(length(U),length(D));
sumU=zeros(size(U));   sigU=zeros(size(U));
sumV=zeros(size(V));   sigV=zeros(size(V));
sumD=zeros(size(D));   sigD=zeros(size(D));

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
    A = exp( 2i*pi*(bsxfun(@times,U*D,wx(freqIds)) + bsxfun(@times,V*D,wy(freqIds))) );
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
    
    
    gradD = sum(bsxfun(@times,gradPx,U) + bsxfun(@times,gradPy,V),1) +.3*gradD;%;%
    gradU = sum(bsxfun(@times,gradPx,D),2) +.3*gradU;%;%
    gradV = sum(bsxfun(@times,gradPy,D),2) +.3*gradV;%;%
    
    gradN2(iter,1) = sum(gradD(:).^2);
    gradN2(iter,2) = sum(gradU(:).^2);
    gradN2(iter,3) = sum(gradV(:).^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Evaluate Residual and stopping criterion          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%residN2(iter) = sum(abs(resid(:)).^2);
%%{
if(mod(iter,EvalPeriod)==0 && ~EvaluatingSigError)
    idEval=iter/EvalPeriod;
    A = exp( 2i*pi*(bsxfun(@times,U*D,wx(freqEvalIds)) + bsxfun(@times,V*D,wy(freqEvalIds))) );
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
        U = sumU/numIterSigError; V = sumV/numIterSigError;D = sumD/numIterSigError;
        sigPx = sqrt(max(0,(sigPx-sumPx.*sumPx/numIterSigError))/numIterSigError);
        sigPy = sqrt(max(0,(sigPy-sumPy.*sumPy/numIterSigError))/numIterSigError);
        sigU = sqrt(max(0,(sigU-sumU.*sumU/numIterSigError))/numIterSigError);
        sigV = sqrt(max(0,(sigV-sumV.*sumV/numIterSigError))/numIterSigError);
        sigD = sqrt(max(0,(sigD-sumD.*sumD/numIterSigError))/numIterSigError);
        break;
    else
        disp(['No convergence after the maximum number of iterations (' num2str(numIterMax) ').']);
        numIterMax = iter + numIterSigError;
        if(numIterSigError>0), EvaluatingSigError = true; end
        disp(['Estimating mean and standard deviation of the parameters (' num2str(numIterSigError) ' more iterations) ...']);
    end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%
%     Line search     %
%%%%%%%%%%%%%%%%%%%%%%%
    stepD = .2/( EpsgradD + sqrt(gradN2(iter,1)));
    stepU = .2/( EpsgradUV + sqrt(sum(gradN2(iter,2:3))));
    stepV = stepU;
    
%%%%%%%%%%%%%%%%%%%%%%%
%       Update        %
%%%%%%%%%%%%%%%%%%%%%%%
    
    D = D - stepD*gradD;
    U = U - stepU*gradU;
    V = V - stepV*gradV;
    
    %%{
    %Normalize U,V and D => the normalized solution is strictly equivlalent w.r.t.
    %the objective function, but the normalization improves the convergence.
    normUV = max(max(U)-min(U),max(V)-min(V))/2;
    D = sort(D) * normUV;
    U = U / normUV;
    V = V / normUV;
    %}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute standard deviation and average parameters with additional iterations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(EvaluatingSigError)
        Px = U*D; Py = V*D;
        sumPx = sumPx + Px; sigPx = sigPx + Px.*Px;
        sumPy = sumPy + Py; sigPy = sigPy + Py.*Py; 
        sumU = sumU + U;    sigU  = sigU + U.*U;
        sumV = sumV + V;    sigV  = sigV + V.*V;
        sumD = sumD + D;    sigD  = sigD + D.*D;
    end


iter=iter+1;
end

Px = gather(Px); sigPx = gather(sigPx);
Py = gather(Py); sigPy = gather(sigPy);
U = gather(U);   sigU = gather(sigU);
V = gather(V);   sigV = gather(sigV);
D = gather(D);   sigD = gather(sigD);

gradN2 = gather(gradN2);
residN2 = gather(residN2);

disp('done!');
end
