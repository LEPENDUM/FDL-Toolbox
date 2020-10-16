%Write YUV image to file with scaling and offset.
%Inputs:
% - filename: Name of the YUV file to read.
% - resX, resY: Horizontal and vertical resolution of the image (not including the padding ensuring image size is a multiple of 8).
% - bitDepth: BitDepth of the YUV file.
% - chromaFormat: '444', '420', '422', or '400'.
% - minVal : Minimum value of the image. This information is used to offset the image (for decoding negative values).
% - nStepQPAdjust: Scaling of pixel values represented by the approximately equivalent adjustment of the HEVC QP.
%
%Outputs:
% - YUV: YUV image.
%

function YUV=readYUV_ScaleForHEVC(filename, resX, resY, bitDepth, chromaFormat, minVal, nStepQPAdjust)

quantizFact = 2^bitDepth-1;
QPQuantLogStep = log(2)/6; %-> an increase of 6 QP steps corresponds a multiplication of the quantization step by 2 within HEVC.
minCUSize = 8;

minVal = double(minVal) / quantizFact;
nStepQPAdjust = double(nStepQPAdjust);

resXPad = minCUSize*ceil(resX/minCUSize);
resYPad = minCUSize*ceil(resY/minCUSize);
YUV = readYUV(filename,resXPad,resYPad,chromaFormat,bitDepth);

%Recover true size (crop padded borders).
YUV = YUV(1:resY,1:resX,:);

%Inverse scale and offset.
valRange = exp(QPQuantLogStep*nStepQPAdjust);

YUV = valRange * double(YUV)/quantizFact + minVal;