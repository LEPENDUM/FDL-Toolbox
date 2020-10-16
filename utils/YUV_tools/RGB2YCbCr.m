%  Conversion function : RGB to YCbCr assuming RGB values are given in gamma (or
%  Log) format in the range [0-max_val]. The Rec-709 chromaticities are used.
%  The returned YCbCr values are also in the range [0-max_val].
function [YCbCr] = RGB2YCbCr(RGB, max_val)

YCbCr(:,:,1,:) = 0.2126 * double(RGB(:,:,1,:)) + 0.7152 * double(RGB(:,:,2,:)) + .0722 * double(RGB(:,:,3,:));
YCbCr(:,:,2,:) = (double(RGB(:,:,3,:)) - YCbCr(:,:,1,:))/1.8556 + max_val/2;
YCbCr(:,:,3,:) = (double(RGB(:,:,1,:)) - YCbCr(:,:,1,:))/1.5748 + max_val/2;