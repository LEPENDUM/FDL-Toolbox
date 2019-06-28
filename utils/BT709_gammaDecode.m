%Apply BT709 standard inverse gamma correction to convert sRGB gamma encoded data to linear BT709.
%(negative input values are clipped to 0).
function Iout = BT709_gammaDecode(Iin)

Mask = Iin < 0.081;
Iout = max(0,Iin./4.5) .* Mask + realpow((max(Iin,0)+0.099)/1.099, 1/0.45) .* (~Mask);