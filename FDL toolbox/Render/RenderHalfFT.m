%GPU Light field Rendering function from FDL representation in the Fourier domain.
%Only half of the spectrum is computed. The other half can be found using symmetries of the Fourier transform assuming a real signal.
%
% Inputs :
% - FDL_gpu : FDL model stored on the GPU. Only the left half of the spectrum is needed. (dimensions : [resY, ceil((1+resX)/2), #color_channels, #layers]).
% - wx_gpu : Grid of horizontal frequency values of all Fourier coefficients (dimensions : [resY, ceil((1+resX)/2)]). wx_gpu must only contain negative or null values (left  half of the spectrum).
% - wy_gpu : Grid of vertical frequency values of all Fourier coefficients (dimensions : [resY, ceil((1+resX)/2])).
% - Disps_gpu : List of disparity values stored on the GPU (dimensions : [1, 1, 1, #layers]).
% - Apfft_gpu : 2D Aperture shape transformed in the disctrete Fourier domain (may only contain left half of the spectrum)-> see buildAperture.
% - dWu, dWv : Horizontal and vertical frequency steps between frequency components in Apfft_gpu -> see buildAperture.
% - uC, vC : Horizontal and vertical indices of the zero frequency component in Apfft_gpu.
% - s : focus value (i.e. disparity value corresponding to the depth in focus).
% - radius : aperture radius.
% - u0, v0 : viewpoint coordinates.

function Ifft = RenderHalfFT(FDL_gpu, wx_gpu, wy_gpu, Disps_gpu, Apfft_gpu, dWu, dWv, uC, vC, s, radius, u0, v0)

%Ifft = sum(arrayfun(@renderLayersCustomAperture,FDL_gpu,wx_gpu,wy_gpu,(s-Disps_gpu)*radius, Disps_gpu),4);
absFact_k = (s-Disps_gpu)*radius;
sgnFact_k = sign(absFact_k);
absFact_k = abs(absFact_k);
Ifft = sum( bsxfun(@times, FDL_gpu, arrayfun(@renderLayersWeightCustomAperture,wx_gpu,wy_gpu, absFact_k, sgnFact_k, Disps_gpu)),4);

    function out = renderLayersWeightCustomAperture(wx,wy, absFact_k, sgnFact_k,disp_k)
%%{
        %Use the absolute value of fact_k=(s-Disp)*radius in order to use
        %the symmetry of the Fourier Transform of real signals ( assuming
        %wx<=0 and dWu>0 we always have u_pos <= uC, so only the left half
        %of the aperture shape spectrum is needed).
        %The conjugate of the result must also be taken when fact_k < 0
        %(i.e. sgnFact_k==-1), in order to exploit the symmetry.
        u_pos = uC + wx*absFact_k/dWu;
        v_pos = vC + wy*absFact_k/dWv;

        if(u_pos>=1 && v_pos>=1 && v_pos<=2*vC-1)
            u_id0 = floor(u_pos);
            v_id0 = floor(v_pos);
            u_id1 = ceil(u_pos);
            v_id1 = ceil(v_pos);
            ur = u_pos-u_id0;
            vr = v_pos-v_id0;

            %Linear interpolation of the aperture function in the Fourier Domain.
            out = Apfft_gpu(v_id0,u_id0)*(1-ur)*(1-vr) + Apfft_gpu(v_id1,u_id0)*(1-ur)*vr + Apfft_gpu(v_id0,u_id1)*ur*(1-vr) + Apfft_gpu(v_id1,u_id1)*ur*vr;
            
            %Take the conjugate when fact_k<0 to exploit symmetry in the Fourier Domain.
            %if(sgnFact_k==-1), out = conj(out);end
            out = real(out) + 1i* sgnFact_k * imag(out);
            
            %Set view position to (u0,v0).
            out = out*exp(2i*pi*disp_k*(u0*wx+v0*wy));
            
        else
            out = single(complex(0));
        end
%}
        %out = single(complex(0));
        %out = exp(2*pi*1i*disp_k*(u0*wx+v0*wy));

    end

end