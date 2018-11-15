%Compute Fourier Disparity Layers using GPU.
%The method uses second order view regularization.
%
%Inputs:
% - ViewsFFT : Fourier Transform of input views with dimensions : (#views, #color channels, Y resolution, X resolution).
% - wx,wy : horizontal and vertical frequency values associated to each frequency component in ViewsFFT. Dimensions = (Y resolution, X resolution).
% - U,V : Vectors of horizontal (U) and vertical (V) view positions.
% - D :  Vector of disparity values.
% - lambda : Regularization parameter.
% - Px : (Optional) Matrix of parameters for the horizontal dimension (default = outer product of U and D).
% - Py : (Optional) Matrix of parameters for the vertical dimension (default = outer product of V and D).
%
%Output:
% FDL model (4D complex array : 2 spatial dimensions, color channels dimension, layers dimension).
% The resulting layers are represented in the Fourier domain.

function FDL = ComputeFDL_gpu(ViewsFFT, wx, wy, U, V, D, lambda,Px,Py)

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
wx = permute(gpuArray(wx(:)),[3 2 1]);
wy = single(wy(:,1:xC));
wy = permute(gpuArray(wy(:)),[3 2 1]);

U = single(U(:));
V = single(V(:));
D = single(D(:)');


if(~exist('Px','var') || ~exist('Py','var') || isempty(Px) || isempty(Py))
    Px = gpuArray(single(U*D));
    Py = gpuArray(single(V*D));
else
    Px = gpuArray(single(Px));
    Py = gpuArray(single(Py));
end

FDL = gpuArray(single(zeros(numDisp,nChan,numFreqProc)));

freqBlockSz = 4096;
%normUV = max(max(U)-min(U),max(V)-min(V))/2;
D = D * 2/(max(D)-min(D));%normalize disparities to adjust the amount of regularization to the disparity range.
d0 = mean(D);

%Adjustment of the regularization parameter to keep the results consistent
%with repect to number of input views and number of layers.
lambdaL2 = .1*numInitViews; %small regularization parameter for l2 norm of the layers.
lambdaViewOrd2Reg = lambda*50*numInitViews*numDisp; %regularization parameter for the 2nd Order view regularization.

for freqSt=1:freqBlockSz:numFreqProc
    freqEnd = min(freqSt+freqBlockSz-1,numFreqProc);
    
    ViewsFFT_gpu = gpuArray(ViewsFFT(:,:,freqSt:freqEnd));
    
    A = exp( 2i*pi*(bsxfun(@times,Px,wx(freqSt:freqEnd)) + bsxfun(@times,Py,wy(freqSt:freqEnd))) );
    AT_ = pagefun(@ctranspose,A);
    
    %%{
    %2nd Order view regularization (infinite camera plane).
    Reg = bsxfun(@times, gpuArray(single(diag((D-d0).^4))), wx(freqSt:freqEnd).^4+wy(freqSt:freqEnd).^4+2*(wx(freqSt:freqEnd).*wy(freqSt:freqEnd)).^2); % 8*pi^4 * ...
    Reg = bsxfun(@plus,lambdaViewOrd2Reg*Reg, lambdaL2*eye(numDisp));
    %}
    
    %{
    %2nd Order view regularization(square portion of the camera plane containing the input viewpoints).
    %Note: if the input viewpoints are not centered on 0, the square portion may not contain them.
    Reg = gpuArray(zeros( numDisp, numDisp, freqEnd-freqSt+1, 'single'));
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
    
    %%{
    %Solve Ridge Regression (linear least squares with Tikhonov Regularization)
    FDL(:,:,freqSt:freqEnd) = ...
        pagefun(@mldivide, ...
            bsxfun(@plus, pagefun(@mtimes,AT_,A), Reg),...
            pagefun(@mtimes,AT_,ViewsFFT_gpu));
    %}
end

clear A AT_ ViewsFFT_gpu



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Find other frequencies with symmetry                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDL(:,:,numFreqProc+1:numFreqY*numFreqX)=0;
FDL = reshape(FDL,numDisp,nChan,numFreqY,numFreqX);
FDL(:,:,yC+1:end,xC:numFreqX) = conj(FDL(:,:,yC-1:-1:1+even_fft(1),xC:-1:1+even_fft(2)));
FDL(:,:,1+even_fft(1):yC,xC+1:numFreqX) = conj(FDL(:,:,end:-1:yC,xC-1:-1:1+even_fft(2)));
FDL = permute(FDL,[3 4 2 1]);

FDL = gather(FDL);

end