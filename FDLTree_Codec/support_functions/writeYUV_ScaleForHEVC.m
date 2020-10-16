%Write YUV image to file with scaling and offset.
%Inputs:
% - filename: Name of the YUV file to write.
% - YUV: Image in YCbCr colorspace but not scaled to the target bitdepth (i.e. floating point format).
% - bitDepth: BitDepth of the YUV file.
% - chromaFormat: '444', '420', '422', or '400'.
% - minVal : Minimum value of the image (after multiplication by 2^bitDepth-1). This information is used to offset the image (for encoding negative values).
% - nStepQPAdjust: Scaling of pixel values represented by the approximately equivalent adjustment of the HEVC QP.
% - append : (default=false) Set to true to append data to the file if it  already exists.
%
%Outputs:
% - resXPad, resYPad: Horizontal and vertical resolution of the yuv file.
% It may be larger than the input image because right and bottom borders
% are padded to ensure the image size is a multiple of 8 (required for
% encoding).
%

function [resXPad, resYPad] = writeYUV_ScaleForHEVC(filename, YUV, bitDepth, chromaFormat, minVal, nStepQPAdjust, append)
if(~exist('append','var'))
    append=false;
end

quantizFact = 2^bitDepth-1;
QPQuantLogStep = log(2)/6; %-> an increase of 6 QP steps corresponds a multiplication of the quantization step by 2 within HEVC.
minCUSize = 8;

minVal = double(minVal) / quantizFact;
nStepQPAdjust = double(nStepQPAdjust);

%Scale and offset according to minVal and nStepQPAdjust parameters:
%minVal is mapped to 0 and a scaling is applied such that encoding the
%scaled version with QP = QP0-nStepQPAdjust is similar to encoding the
%unscaled version encoded with QP=QP0 (the difference is the internal precision for the same bitDetph used).
valRange = exp(QPQuantLogStep*nStepQPAdjust);
YUV = round( quantizFact * (double(YUV) - minVal) / valRange );
if( bitDepth < 1 || bitDepth > 16 )
    error('Unsupported bitdepth: must be between 1 and 16.');
elseif(bitDepth <= 8)
    YUV = uint8(YUV);
elseif(bitDepth <= 16)
    YUV = uint16(YUV);
end

%Pad borders to enforce a size multiple of the minimum CU Size of HEVC.
resX = size(YUV,2);
resY = size(YUV,1);
resXPad = minCUSize*ceil(resX/minCUSize);
resYPad = minCUSize*ceil(resY/minCUSize);
YUV(:,end+1:end+resXPad-resX,:) = 2^(bitDepth-1);
YUV(end+1:end+resYPad-resY,:,:) = 2^(bitDepth-1);

writeYUV(YUV,filename,chromaFormat,append);