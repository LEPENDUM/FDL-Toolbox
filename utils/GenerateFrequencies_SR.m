% Generate array of frequency values for super-resolution of an image of
% input resolution 'resX', 'resY' and an integer scaling factor 'scale'.
% The frequencies are split into scaleX*scaleY regions, where the regions
% are defined to perform the convolution needed for super-resolution.
%
% Inputs:
% -resX,resY: Input sampling resolution.
% -scale:     Scaling factor(s) for super-resolution (must be a strictly positive integer).
%             Scale can be a scalar or a vector of 2x1 containing horizontal and vertiacal scaling factors.
% -hexaSampling (Optional): 
%       0 (default) -> Square input sampling.
%       1->Hexagonal input sampling with odd horizontal lines shifted to the right (unknown top left corner pixel).
%       2->Hexagonal input sampling with even horizontal lines shifted to the right (known top left corner pixel).
%       For Hexagonal samplings, the vertical resolution must be a multiple of 2
%        and the horizontal scale must be either 1 or a multiple of 2.
% - halfSpectrum (Optional):
%       true (default) -> Generate only left half of the frequency spectrum
% (only half of the frequencies are needed to process real signals due to symmetries).
%       false-> Generate full frequency spectrum.
%
% Outputs:
% -wx : Array of horizontal frequency values.
% -wy : Array of vertical frequency values.
% -xC : Horizontal index of the zero frequency ( wx(:,xC,c)=0 with c the index of the central region / xAxis oriented to the right).
% -yC : Vertical index of the zero frequency ( wy(yC,:,c)=0 with c the index of the central region / yAxis oriented upwards).
% -xRanges :   Horizontal indices of the frequencies of each region in the (half) frequency spectrum.
% -yRanges :   Vertical indices of the frequencies of each region in the frequency spectrum.
% -coeffs :    Coefficient associated to each region for the Fourier domain convolution (spatial domain masking).
% -centralFreqsId : Index of the region correponding to low frequencies.
% -regionIds : (Optional) Indices of the frequencies of each region in the (half) frequency spectrum.

function [wx,wy,xC,yC,xRanges,yRanges,coeffs,centralFreqsId,regionIds] = GenerateFrequencies_SR(resX,resY,scale,hexaSampling,halfSpectrum)

    scaleX = scale(1);
    if(numel(scale)==1)
        scaleY = scaleX;
    else
        scaleY = scale(2);
    end

    if(~exist('hexaSampling','var')),hexaSampling=0;end
    if(~exist('halfSpectrum','var')),halfSpectrum=true;end
    
    
    if(any(round(scale)-scale~=0) || any(scale<=0) )
        error('scale factor must be a strictly positive integer value.');
    end
    if(hexaSampling>0 && ((mod(scaleX,2)~=0 && scaleX~=1) || mod(resY,2)~=0) )
%        error('Vertical resolution ''resY'' and horizontal scale ''scaleX'' must be multiples of 2 for hexagonal samplings.');
    end
    
    resXs = scaleX*resX;
    resYs = scaleY*resY;
    
    centralFreqsId = sub2ind([scaleY scaleX],ceil((scaleY+1)/2),ceil((scaleX+1)/2));
    numCoeffs = scaleX*scaleY;
    
    if(hexaSampling==0) 
        %Square sampling
        [convShiftsX,convShiftsY] = meshgrid(1:scaleX,1:scaleY);
        convShiftsX = (convShiftsX(:)-ceil((scaleX+1)/2)) * resX;
        convShiftsY = (convShiftsY(:)-ceil((scaleY+1)/2)) * resY;
        coeffs = ones(numCoeffs,1);
    else
        %Hexagonal sampling
        %{
        %[convShiftsY2,convShiftsX2]=ind2sub([2*stepY,stepX],2-mod(stepY,2):2:2*stepY*stepX);
        %convShiftsX2 = (convShiftsX2(:)-ceil((stepX+1)/2)) * imgSizeIn(2);
        %convShiftsY2 = (convShiftsY2(:)-(stepY+1)) * floor(imgSizeIn(1)/2);
        [convShiftsY2,convShiftsX2]=ind2sub([scaleY*2,scaleX],find(M_mini_f~=0));
        convShiftsX2 = (convShiftsX2(:)-ceil((scaleX+1)/2)) * imgSizeIn(2);
        convShiftsY2 = (convShiftsY2(:)-(scaleY+1)) * floor(imgSizeIn(1)/2);
        %}
        M_mini = zeros(scaleY*2,scaleX);
        M_mini(1:2*scaleY:end,1+floor(scaleX/2)*(2-hexaSampling):scaleX:end)=1;
        M_mini(1+scaleY:2*scaleY:end,1+floor(scaleX/2)*(hexaSampling-1):scaleX:end)=1;
        M_mini_f = fftshift(fft2(M_mini));
        coeffs = M_mini_f(M_mini_f~=0)/2;
        [convShiftsY,convShiftsX]=ind2sub([scaleY*2,scaleX],find(M_mini_f~=0));
        convShiftsX = (convShiftsX(:)-ceil((scaleX+1)/2)) * resX;
        convShiftsY = (convShiftsY(:)-(scaleY+1)) * resY/2;
    end
    
    [wx0,wy0] = meshgrid(1:resXs,1:resYs);
    xCs = ceil((resXs+1)/2);
    yCs = ceil((resYs+1)/2);
    wx0=single((wx0-xCs)/resXs);
    wy0=single((yCs-wy0)/resYs);
    xC = ceil((resX+1)/2); 
    yC = ceil((resY+1)/2);
    
    yC_b = floor((resY-1)/2);
    
    if(halfSpectrum)
        %Generate only half of the frequencies (only one half spectrum necessary to process real signals).
        wx = zeros(resY,xC,numCoeffs,'single');
        wy = zeros(resY,xC,numCoeffs,'single');
        xRanges = zeros(xC,numCoeffs);
        yRanges = zeros(resY,numCoeffs);
        range_x_center = 1+xCs-xC : xCs;
        range_y_center = 1+yCs-yC : yCs+yC_b;
    else
        %Generate all the frequencies
        xC_r = floor((resX-1)/2);
        wx = zeros(resY,resX,numCoeffs,'single');
        wy = zeros(resY,resX,numCoeffs,'single');
        xRanges = zeros(resX,numCoeffs);
        yRanges = zeros(resY,numCoeffs);
        range_x_center = 1+xCs-xC : xCs+xC_r;
        range_y_center = 1+yCs-yC : yCs+yC_b;
    end
    
    for i=1:numCoeffs
        xRanges(:,i) = mod( range_x_center + convShiftsX(i) - 1, resXs)+1;
        yRanges(:,i) = mod( range_y_center + convShiftsY(i) - 1, resYs)+1;
        
        wx(:,:,i) = wx0(yRanges(:,i),xRanges(:,i));
        wy(:,:,i) = wy0(yRanges(:,i),xRanges(:,i));
    end
   
    if(nargout>8)
        regionIds = zeros(numCoeffs,size(xRanges,1)*resY);
        for i=1:numCoeffs
            [idX,idY]=meshgrid(xRanges(:,i),yRanges(:,i));
            regionIds(i,:)=sub2ind([resYs resXs],idY(:),idX(:));
        end
    end
    
end