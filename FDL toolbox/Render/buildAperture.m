%Function generating an aperture shape image and its Fourier domain representation.
%
%Inputs :
% - shape : (string) 'disk', 'ring', 'rect', 'polygon', 'dirac'
% - rad : radius in number of pixels. Larger radius increases precision in the pixel domain.
% - pad : padding size before computing the Fourier transform. Padding increases precision in the discrete Fourier domain.
% - extraParams : additional parameters specific to some of the shapes.
%     -#1: thickness : thickness between 0 (i.e. thin border) and 1 (entierly fills the aperture area). Only implemented for ring and polygon shapes.
%     -#2: numBlades : number of segments for the polygon shape.
%     -#3: startAngle: rotation angle for the polygon shape.

%Outputs :
% - Ap : Aperture shape in the pixel domain.
% - ApUnitFFTHalf : Aperture shape transformed in the disctrete Fourier domain. Only returns the left half of the spectrum (the rest is obtained by symmetry of the complex conjugate).
% - UnitDwx, UnitDwy : Horizontal and vertical frequency steps between frequency components in ApUnitFFTHalf assuming the generated aperture shape has a unit radius.
% - TrCX, TrCY : Horizontal and vertical indices of the zero frequency component in ApUnitFFTHalf.
% - TrueRadius : Corrected value of the input radius after drawing the shape Ap.
% The true radius compensates for possible inaccuracies or inconsistencies in the
% drawing functions.
%
% Note #1: The Fourier domain representation (ApUnitFFTHalf, UnitDwx, UnitDwy, TrCX,TrCY)
% results in an aperture with normalized radius. Hence the input radius
% parameter does not affect the scale, but is only used to draw an
% unaliased shape with unit radius.
%
% Note #2: In Matlab's convention for image indexing, the vertical axis is
% directed downwards. This function returns a negative value for UnitDwy to use the
% convention with vertical axis directed upwards.

function [Ap,ApUnitFFTHalf,UnitDwx,UnitDwy,TrCX,TrCY,TrueRadius] = buildAperture(shape,rad,pad,extraParams)

szApX = 2*rad+1;
szApY = 2*rad+1;
szTrX = 2*(pad+rad)+1;
szTrY = 2*(pad+rad)+1;
TrCX = pad+rad+1;
TrCY = pad+rad+1;

if(exist('extraParams','var') && ~isempty(extraParams))
    thickness = extraParams(1);
    extraParams=extraParams(2:end);
else
    extraParams=[];
    thickness = 1;
end


switch shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'disk'
Ap = insertShape(zeros(szApY,szApX),'FilledCircle',[rad+1 rad+1 rad],'Color','white','LineWidth',1,'Opacity',1);
%A(size+2:end,1:size)=A(size:-1:1,1:size);A(1:size,size+2:end)=A(1:size,size:-1:1);A(size+2:end,size+2:end)=A(size:-1:1,size:-1:1);
Ap = Ap(:,:,1);
Apad = padarray(Ap,[pad pad]);
ApUnitFFTHalf = fftshift(fft2(ifftshift(Apad)));
TrueRadius = sqrt(abs(ApUnitFFTHalf(TrCY,TrCX)) / pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'ring'
lineThickness = max(min(round(1+rad*thickness),rad+1),2);
Ap = insertShape(zeros(szApY,szApX),'Circle',[rad+1 rad+1 rad-floor(lineThickness/2)],'Color','white','LineWidth',lineThickness);
TrueRadius = rad+.51;
Ap = Ap(:,:,1);
Apad = padarray(Ap,[pad pad]);
ApUnitFFTHalf = fftshift(fft2(ifftshift(Apad)));
%TrueThickness = TrueRadius - sqrt(TrueRadius*TrueRadius - abs(ApUnitFFTHalf(TrCY,TrCX))/pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'rect'
%By convention, the radius is half the rectangle size.
Ap = insertShape(zeros(szApY,szApX),'FilledRectangle',[1 1 2*rad+1 2*rad+1],'Color','white','LineWidth',1,'Opacity',1);
TrueRadius = rad+.5;
Ap = Ap(:,:,1);
Apad = padarray(Ap,[pad pad]);
ApUnitFFTHalf = fftshift(fft2(ifftshift(Apad)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'polygon'
%By convention, the polygon is regular and inscribed in a circle of the input radius.

if(length(extraParams)>0), numBlades = extraParams(1); else, numBlades=6;end
if(length(extraParams)>1), startAngle = extraParams(2); else, startAngle=0;end
if(numBlades<2), error('Polygon aperture cannot have less than 2 vertices.');end

angles = [0:2*pi/numBlades:2*pi*(numBlades-1)/numBlades];
angles = angles + pi/2-pi/numBlades + startAngle;
pos = zeros(1,2*numBlades);
pos(1:2:end) = rad+1 + rad*cos(angles);%X position of polygon vertices.
pos(2:2:end) = rad+1 + rad*sin(angles);%Y position of polygon vertices.
 
Ap = insertShape(zeros(szApY,szApX),'FilledPolygon',pos,'Color','white','LineWidth',1,'Opacity',1);
Ap = Ap(:,:,1);

F=sqrt((1-cos(2*pi/numBlades))/2);
F=numBlades*F*sqrt(1-F^2);
TrueRadius = sqrt(sum(Ap(:))/F);

if(thickness<1)
    radInterior = min(rad-2,max(0, (rad-1) * (1-thickness)));
    pos(1:2:end) = rad+1 + radInterior*cos(angles);
    pos(2:2:end) = rad+1 + radInterior*sin(angles);
    ApInterior = insertShape(zeros(szApY,szApX),'FilledPolygon',pos,'Color','white','LineWidth',1,'Opacity',1);
    ApInterior = ApInterior(:,:,1);
    Ap = Ap - 1*ApInterior;
end

Apad = padarray(Ap,[pad pad]);
ApUnitFFTHalf = fftshift(fft2(ifftshift(Apad)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'dirac'
Ap = zeros(szApY,szApX);Ap(rad+1,rad+1)=1;
ApUnitFFTHalf = ones(szTrY,szTrX);
TrueRadius = 1;
    otherwise
Ap = zeros(szApY,szApX);Ap(rad+1,rad+1)=1;
ApUnitFFTHalf = ones(szTrY,szTrX);
TrueRadius = 1;
end

%ApUnitFFTHalf = ApUnitFFTHalf/ApUnitFFTHalf(TrCY,TrCX);
ApUnitFFTHalf = ApUnitFFTHalf(:,1:TrCX);
UnitDwx = TrueRadius/(szTrX-1);
UnitDwy = -TrueRadius/(szTrY-1); %Invert the sign to use the convention of a vertical axis pointing upwards ('Bottom->Up' convention).