%CPU Light field Rendering function from FDL representation in the Fourier domain.
%Only half of the spectrum is computed. The other half can be found using symmetries of the Fourier transform assuming a real signal.
%
% Inputs :
% - FDL : FDL model. Only the left half of the spectrum is needed. (dimensions : [resY, ceil((1+resX)/2), #color_channels, #layers]).
% - wx : Grid of horizontal frequency values of all Fourier coefficients (dimensions : [resY, ceil((1+resX)/2)]). wx must only contain negative or null values (left  half of the spectrum).
% - wy : Grid of vertical frequency values of all Fourier coefficients (dimensions : [resY, ceil((1+resX)/2])).
% - Disps : List of disparity values (dimensions : [1, 1, 1, #layers]).
% - Apfft : 2D Aperture shape transformed in the disctrete Fourier domain (may only contain left half of the spectrum)-> see buildAperture.
% - dWu, dWv : Horizontal and vertical frequency steps between frequency components in Apfft -> see buildAperture.
% - uC, vC : Horizontal and vertical indices of the zero frequency component in Apfft.
% - s : focus value (i.e. disparity value corresponding to the depth in focus).
% - radius : aperture radius.
% - u0, v0 : viewpoint coordinates.

function Ifft = RenderHalfFT_cpu(FDL, wx, wy, Disps, Apfft, dWu, dWv, uC, vC, s, radius, u0, v0)

absFact = (s-Disps)*radius;
sgnFact = sign(absFact);
absFact = abs(absFact);

resY = size(FDL,1);
resX = size(FDL,2);
numDisp = size(FDL,4);
Weigths = zeros(resY,resX,1,numDisp);

for k = 1:numDisp
    %Use the absolute value of fact_k=(s-Disp)*radius in order to use
    %the symmetry of the Fourier Transform of real signals ( assuming
    %wx<=0 and dWu>0 we always have u_pos <= uC, so only the left half
    %of the aperture shape spectrum is needed).
    %The conjugate of the result must also be taken when fact_k < 0
    %(i.e. sgnFact_k==-1), in order to exploit the symmetry.
    u_pos = uC + wx*absFact(k)/dWu;
    v_pos = vC + wy*absFact(k)/dWv;
    
    %Linear interpolation of the aperture function in the Fourier Domain.
    Weigths(:,:,1,k) = interp2(Apfft,u_pos,v_pos,'linear',0);
    
    %Take the conjugate when fact_k<0 to exploit symmetry in the Fourier Domain.
    if(sgnFact(k)==-1), Weigths(:,:,1,k) = conj(Weigths(:,:,1,k));end
end

%Set view position to (u0,v0).
Weigths = bsxfun(@times, Weigths, exp(2i*pi*bsxfun(@times,Disps, u0*wx+v0*wy)) );

%Apply weights (shift & filter) and sum the layers
Ifft = sum( bsxfun(@times, FDL, Weigths),4);

end