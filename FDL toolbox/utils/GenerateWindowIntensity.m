%Generate intensity profile for windowing.
% Inputs:
% - resX : full horizontal image resolution.
% - resY : full vertical image resolution.
% - borderX : Number of pixels with intensity decay on the left and right borders.
% - borderY : Number of pixels with intensity decay on the top and bottom borders.
% - windowProfile : 'hann' or 'linear' window profiles.
%
% Output:
% -Window : Window intensity profile. Only the pixels on the borders have intensity inferior to 1.

function [Window] = GenerateWindowIntensity(resX,resY,borderX,borderY, windowProfile)
    [gridx,gridy] = meshgrid(1:resX,1:resY);
    borderXSup1 = max(1,borderX);
    borderYSup1 = max(1,borderY);
    
    xW = borderXSup1*ones(1,resX);
    yW = borderYSup1*ones(1,resY);
    
    
    if(strcmp(windowProfile,'linear'))
        xW(1:borderX) = 0:borderX-1;
        xW(end+1-borderX:end) = borderX-1:-1:0;
        xW = xW / borderXSup1;
        
        yW(1:borderY) = 0:borderY-1;
        yW(end+1-borderY:end) = borderY-1:-1:0;
        yW = yW / borderYSup1;
    
    elseif(strcmp(windowProfile,'hann'))
        xW(1:borderX) = 0:borderX-1;
        xW(end+1-borderX:end) = 1+borderX:2*borderX;
        xW = (.5-.5*cos(pi*xW/borderXSup1)).*(xW>=0 & xW<=2*borderXSup1);
                
        yW(1:borderY) = 0:borderY-1;
        yW(end+1-borderY:end) = 1+borderY:2*borderY;
        yW = (.5-.5*cos(pi*yW/borderYSup1)).*(yW>=0 & yW<=2*borderYSup1);
    else
        error('Unknown window profile');
    end
    
    Window = single(xW(gridx).*yW(gridy));
    

end