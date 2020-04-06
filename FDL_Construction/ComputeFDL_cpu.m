%Compute Fourier Disparity Layers using CPU.
%The method uses second order view regularization.
%
%Inputs:
% - ViewsFFT : Fourier Transform of input views as a 4D array with dimensions : 1.Views, 2.Color channels, 3.Vertical axis (Y), 4.Horizontal axis (X).
% - wx,wy : horizontal and vertical frequency values associated to each frequency component in ViewsFFT. size = (Y resolution, X resolution).
% - U,V : Vectors of horizontal (U) and vertical (V) view positions.
% - D :  Vector of disparity values.
% - lambdas : (Optional) Regularization parameters. Can be either empty, scalar or a vector of 2 elements:
%       -lambdas(1): Regularization parameter controling l2 regularization term (1 by default).
%       -lambdas(2): Regularization parameter controling 2nd order view regularization (same as lambdas(1) by default).
% - Px : (Optional) Matrix of parameters for the horizontal dimension (default = outer product of U and D).
% - Py : (Optional) Matrix of parameters for the vertical dimension (default = outer product of V and D).
% - varargin: Additional optional arguments given as (Name,Value) pairs:
%   - 'freqBlockSz' (default=4096) : Number of frequencies processed in parallel.
%   - 'ViewWeights' (default=[]): Weights of each view for the linear regression (lower weight=lower confidence in a view) => vector (#views x 1).
%   If it is empty no weights are used (equivalent to all the weights equal to 1).
%
%Output:
% FDL model (4D complex array with dimensions: 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color channels, 4.Layers.
% The resulting layers are represented in the Fourier domain.

function FDL = ComputeFDL_cpu(ViewsFFT, wx, wy, U, V, D, lambdas, Px, Py, varargin)

p = inputParser;
addParameter (p,'ViewWeights',  [],   @isnumeric);
addParameter (p,'freqBlockSz',  4096, @isnumeric);
parse(p,varargin{:});
ViewWeights = p.Results.ViewWeights;
freqBlockSz = p.Results.freqBlockSz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nChan = size(ViewsFFT,2);
xC = ceil((size(ViewsFFT,4)+1)/2);
yC = ceil((size(ViewsFFT,3)+1)/2);
even_fft = 1-mod([size(ViewsFFT,3) size(ViewsFFT,4)],2);

numInitViews = size(ViewsFFT,1);
numDisp = length(D);
numFreqX = size(ViewsFFT,4);
numFreqY = size(ViewsFFT,3);
numFreqProc = numFreqY*(xC-1)+yC;

ViewsFFT = single(ViewsFFT);
wx = single(wx(:,1:xC));
wx = permute(wx(:),[3 2 1]);
wy = single(wy(:,1:xC));
wy = permute(wy(:),[3 2 1]);

U = single(U(:));
V = single(V(:));
D = single(D(:)');

if(~exist('Px','var') || ~exist('Py','var') || isempty(Px) || isempty(Py))
    Px = single(U*D);
    Py = single(V*D);
else
    Px = single(Px);
    Py = single(Py);
end

%View Weights parameters (weight the different views in the FDL construction).
useViewWeights = ~isempty(ViewWeights);
if(useViewWeights)
    if(numel(ViewWeights)~=numInitViews)
        error('If ''ViewWeights'' is not empty, the number of elements should be equal to the number of views.');
    end
    ViewWeights = reshape(single(ViewWeights),1,[]).^2;
    useViewWeights=true;
end

if(~exist('lambdas','var') || isempty(lambdas))
    lambdas = [1, 1];
elseif(numel(lambdas)==1)
    lambdas(2)=lambdas(1);
end

FDL = single(zeros(numDisp,nChan,numFreqProc));

%normUV = max(max(U)-min(U),max(V)-min(V))/2;
D = D * 2/(max(D)-min(D));%normalize disparities to adjust the amount of regularization to the disparity range.
d0 = mean(D);

%Adjustment of the regularization parameter to keep the results consistent
%with repect to number of input views and number of layers.
lambdaL2 = lambdas(1)*numInitViews*numDisp*.005; %small regularization parameter for l2 norm of the layers.
lambdaViewOrd2Reg = lambdas(2)*50*numInitViews*numDisp; %regularization parameter for the 2nd Order view regularization.

for freqSt=1:freqBlockSz:numFreqProc
    freqEnd = min(freqSt+freqBlockSz-1,numFreqProc);
    
    A = exp( 2i*pi*(bsxfun(@times,Px,wx(freqSt:freqEnd)) + bsxfun(@times,Py,wy(freqSt:freqEnd))) );
    
    %%{
    %2nd Order view regularization (infinite camera plane).
    Reg = bsxfun(@times, single(diag((D-d0).^4)), wx(freqSt:freqEnd).^4+wy(freqSt:freqEnd).^4+2*(wx(freqSt:freqEnd).*wy(freqSt:freqEnd)).^2); % 8*pi^4 * ...
    Reg = bsxfun(@plus,lambdaViewOrd2Reg*Reg, lambdaL2*eye(numDisp));
    %}
    
    %{
    %2nd Order view regularization(square portion of the camera plane containing the input viewpoints).
    %Note: if the input viewpoints are not centered on 0, the square portion may not contain them.
    Reg = zeros( numDisp, numDisp, freqEnd-freqSt+1, 'single');
    for j=1:numDisp
        for k=1:j
            Reg(j,k,:) = ((D(k)-d0).*(D(j)-d0)).^2 * ...
                        (wx(freqSt:freqEnd).^4+wy(freqSt:freqEnd).^4+2*(wx(freqSt:freqEnd).*wy(freqSt:freqEnd)).^2) .* ...
                        sinc(2*normUV*(D(k)-D(j))*wx(freqSt:freqEnd)) .* sinc(2*normUV*(D(k)-D(j))*wy(freqSt:freqEnd));
            Reg(k,j,:) = Reg(j,k,:);
        end
    end
    Reg = bsxfun(@plus,lambdaViewOrd2Reg*Reg, lambdaL2*eye(numDisp));
    %}
    
    if(useViewWeights)
        for i = 1:freqEnd-freqSt+1
            freqId = freqSt + i - 1;
            ATw = bsxfun(@times, A(:,:,i)',ViewWeights);
            FDL(:,:,freqId) = (ATw * A(:,:,i) + Reg(:,:,i)) \ ( ATw * ViewsFFT(:,:,freqId));
        end
    else
        for i = 1:freqEnd-freqSt+1
            freqId = freqSt + i - 1;
            FDL(:,:,freqId) = (A(:,:,i)'*A(:,:,i) + Reg(:,:,i)) \ ( A(:,:,i)' * ViewsFFT(:,:,freqId));
        end
    end
end

clear A ATw ViewsFFT_gpu



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Find other frequencies with symmetry                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDL(:,:,numFreqProc+1:numFreqY*numFreqX)=0;
FDL = reshape(FDL,numDisp,nChan,numFreqY,numFreqX);
FDL(:,:,yC+1:end,xC:numFreqX) = conj(FDL(:,:,yC-1:-1:1+even_fft(1),xC:-1:1+even_fft(2)));
FDL(:,:,1+even_fft(1):yC,xC+1:numFreqX) = conj(FDL(:,:,end:-1:yC,xC-1:-1:1+even_fft(2)));
FDL = permute(FDL,[3 4 2 1]);

end