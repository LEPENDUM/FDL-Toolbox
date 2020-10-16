%  Conversion function : YCbCr to RGB 
%  The returned RGB values are gamma corrected and in the range [0-max_val]. The Rec-709 chromaticities are used.
function [RGB] = YCbCr2RGB(YCbCr, max_val)

RGB(:,:,1,:) = 1.5748*(double(YCbCr(:,:,3,:))-(max_val/2)) + double(YCbCr(:,:,1,:));
RGB(:,:,3,:) = 1.8556*(double(YCbCr(:,:,2,:))-(max_val/2)) + double(YCbCr(:,:,1,:));
RGB(:,:,2,:) = (double(YCbCr(:,:,1,:)) - 0.2126 * RGB(:,:,1,:) - .0722 * RGB(:,:,3,:)) / 0.7152;