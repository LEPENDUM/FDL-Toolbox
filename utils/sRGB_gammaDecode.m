%Apply sRGB standard inverse gamma correction to convert sRGB gamma encoded data to linear RGB.
%(negative input values are clipped to 0).
function Iout = sRGB_gammaDecode(Iin)

Mask = Iin < 0.04045;
Iout = max(0,Iin/12.92) .* Mask + realpow((max(Iin,0)+0.055)/1.055,2.4) .* (~Mask);