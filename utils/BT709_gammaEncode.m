%Apply BT709 standard gamma correction to convert from linear RGB data to BT709.
%(negative input values are clipped to 0).
function Iout = BT709_gammaEncode(Iin)

Mask = Iin < 0.018;
Iout = max(0,Iin * 4.5) .* Mask + (1.099*realpow(max(Iin,0),0.45) - 0.099) .* (~Mask);