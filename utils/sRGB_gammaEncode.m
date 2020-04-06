%Apply sRGB standard gamma correction to convert from linear RGB data to sRGB.
%(negative input values are clipped to 0).
function Iout = sRGB_gammaEncode(Iin)

Mask = Iin < 0.0031308;
Iout = max(0,Iin*12.92) .* Mask + (1.055*realpow(max(Iin,0),1/2.4) - 0.055) .* (~Mask);