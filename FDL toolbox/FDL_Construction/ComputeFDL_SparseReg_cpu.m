%Compute Fourier Disparity Layers using CPU
%The method uses ElasticNet regularization (both l1 and l2).
%
%Inputs:
% - ViewsFFT : Fourier Transform of input views with dimensions : (#views, #color channels, Y resolution, X resolution).
% - wx,wy : horizontal and vertical frequency values associated to each frequency component in ViewsFFT. Dimensions = (Y resolution, X resolution).
% - U,V : Vectors of horizontal (U) and vertical (V) view positions.
% - D :  Vector of disparity values.
% - lambdaL1 : Regularization parameter for l1 norm penalty term.
% - lambdaL2 : Regularization parameter for l2 norm penalty term.
% - Px : (Optional) Matrix of parameters for the horizontal dimension (default = outer product of U and D).
% - Py : (Optional) Matrix of parameters for the vertical dimension (default = outer product of V and D).
%
%Output:
% FDL model (4D complex array : 2 spatial dimensions, color channels dimension, layers dimension).
% The resulting layers are represented in the Fourier domain.

function FDL = ComputeFDL_SparseReg_cpu(ViewsFFT, wx, wy, U, V, D, lambdaL1, lambdaL2,Px,Py)

nChan = size(ViewsFFT,2);
xC = ceil((size(ViewsFFT,4)+1)/2);
yC = ceil((size(ViewsFFT,3)+1)/2);
even_fft = 1-mod([size(ViewsFFT,3) size(ViewsFFT,4)],2);

numDisp = length(D);
numFreqX = size(ViewsFFT,4);
numFreqY = size(ViewsFFT,3);
numFreqProc = numFreqY*(xC-1)+yC;


ViewsFFT = single(ViewsFFT);
wx = single(wx(:,1:xC));
wy = single(wy(:,1:xC));

U = single(U(:));
V = single(V(:));
D = single(D(:)');


if(~exist('Px','var') && ~exist('Py','var'))
    Px = U*D;
    Py = V*D;
end

FDL = single(zeros(numDisp,nChan,numFreqProc));

if(lambdaL1==0)
    Reg = single(lambdaL2*eye(numDisp));
    parfor freq=1:numFreqProc
        b = ViewsFFT(:,:,freq);
        A = exp( 2i*pi*(wx(freq)*Px + wy(freq)*Py) );
        FDL(:,:,freq) = (A'*A + Reg) \ ( A' * b);
    end
else
    parfor freq=1:numFreqProc
        b = ViewsFFT(:,:,freq);
        A = exp( 2i*pi*(wx(freq)*Px + wy(freq)*Py) );
        [FDL(:,:,freq),~] = elasticNetADMM(b,A,lambdaL1,lambdaL2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Find other frequencies with symmetry                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDL(:,:,numFreqProc+1:numFreqY*numFreqX)=0;
FDL = reshape(FDL,numDisp,nChan,numFreqY,numFreqX);
FDL(:,:,yC+1:end,xC:numFreqX) = conj(FDL(:,:,yC-1:-1:1+even_fft(1),xC:-1:1+even_fft(2)));
FDL(:,:,1+even_fft(1):yC,xC+1:numFreqX) = conj(FDL(:,:,end:-1:yC,xC-1:-1:1+even_fft(2)));
FDL = permute(FDL,[3 4 2 1]);


end