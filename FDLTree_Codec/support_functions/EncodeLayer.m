%Function for encoding a layer.
%
%Inputs:
% - Layer : Layer image to encode.
%
% - writeFileFmt: filename for the output files. Suffixe and extension are added depending on the generated files.
%
% - refLayerFileFmt: Filename of the reference yuv file for encoding the current layer. May be left empty of the layer is intra coded.
%
% - QPVirtual :
%      HEVC encoding QP parameter assuming no rescaling operation applies.
%      If the data is scaled (to cover the range of values in the given bitdepth), the QP is adjusted to obtain similar overall amount of quantization than without the scaling and using the virtual QP.
%      Scaling and QP adjustment are done to increase the internal encoder's precision without changing the bitdepth
%
% - cfg.chromaFormat : YUV Chroma format (default='444'):
%                      - '444' encode for full resolution chroma.
%                      - '420' chroma downsampling.
%                      - '400' only Luma.
% - cfg.bitDepth     : Encoding bitdepth in HEVC.
% - cfg.splitCodec   : Codec used for Split Images:
%                      - 'HEVC' : standard HEVC (only intra).
%                      - 'HEVC+linILP' (default) : HEVC with additional block-wise linear Inter Layer Prediction mode.
%                      - 'HEVC+nonlinILP' : HEVC with additional block-wise non-linear Inter Layer Prediction mode.
% - cfg.skipYUVConv    (default=false) : Set true to skip the YUV conversions in the encoder (e.g. if the input FDL is already given in YUV colorspace).
% - cfg.lFact          (default=1)     : Factor to multiply the lagrangian parameter lambda used in HEVC for rate distortion optimisation (only used if cfg.splitCodec is 'HEVC+linILP' or 'HEVC+nonlinILP').
% - cfg.outputModePred (default=false) : Set true to generate prediction image and mode decision image (only used if cfg.splitCodec is 'HEVC+linILP' or 'HEVC+nonlinILP').
%
%Outputs:
% - LayerRec: Reconstructed layer.
% - nBits: Number of bits of encoded image.
% - outFiles: Structure containing the names of all the output files.
%       - outFiles.log -> log file.
%       - outFiles.org -> yuv file before compression.
%       - outFiles.rec -> yuv file file reconstructed after compression.
%       - outFiles.bin -> binary file (compressed layer).
%       - outFiles.ext -> structure containing suffix and extension for each of the log/org/rec/bin files.
% - minLayVal: Minimum value of the layer (after multiplication by  2^cfg.bitDepth-1). This information is required to read the image with correct offset.
% - QPAdjust: Adjustment of the encoding QP compensated by scaling pixel values. This information is required to read the image with correct scale.
% - MSE: Mean square error of the reconstructed image.
%

function [LayerRec, nBits, outFiles, minLayVal, QPAdjust, MSE] = EncodeLayer(Layer, writeFileFmt, refLayerFileFmt, QPVirtual, cfg)
resX = size(Layer,2);
resY = size(Layer,1);
bitDepth = cfg.bitDepth;
if(isfield(cfg,'chromaFormat')),    chromaFormat=cfg.chromaFormat;       else, chromaFormat='444';end
if(isfield(cfg,'skipYUVConv')),     skipYUVConv=cfg.skipYUVConv;         else, skipYUVConv=false;end
if(isfield(cfg,'splitCodec')),      splitCodec = cfg.splitCodec;         else, splitCodec='HEVC+linILP';end
if(isfield(cfg,'lFact')),           lFact=cfg.lFact;                     else, lFact=1;end
if(isfield(cfg,'outputModePred')),  outputModePred=cfg.outputModePred;   else, outputModePred=false;end


if(~exist('refLayerFileFmt','var') || isempty(refLayerFileFmt))
    useRef = false;
    HEVC_exe = [HEVC_ILP_path 'TAppEncoder_HM-12.1+RExt-4.2.exe'];
else
    refLayerFileFmt = fullfile(refLayerFileFmt);
    refLayerFile = [refLayerFileFmt '_REC.yuv'];
    useRef = true;
    
    switch splitCodec
        case 'HEVC'
            useRef = false;
            HEVC_exe = [HEVC_ILP_path 'TAppEncoder_HM-12.1+RExt-4.2.exe'];
        case 'HEVC+linILP'
            HEVC_exe = [HEVC_ILP_path 'TAppEncoder_ILP_lin_trFast.exe'];
        case 'HEVC+nonlinILP'
            HEVC_exe = [HEVC_ILP_path 'TAppEncoder_ILP_3lin_Ext_HybridFull.exe'];
        otherwise
            error(['Unknown encoding type ''' splitCodec ''' specified in the splitCodec configuration field.']);
    end
end

HEVC_cfg = [HEVC_ILP_path 'encoder_intra_main10.cfg'];

writeFileFmt = fullfile(writeFileFmt);

outFiles.ext.log = '.log';
outFiles.ext.org = '_ORG.yuv';
outFiles.ext.rec = '_REC.yuv';
outFiles.ext.bin = '.bin';
outFiles.log = [writeFileFmt outFiles.ext.log];
outFiles.org = [writeFileFmt outFiles.ext.org];
outFiles.rec = [writeFileFmt outFiles.ext.rec];
outFiles.bin = [writeFileFmt outFiles.ext.bin];

if(outputModePred)
    outFiles.ext.mode = '_MODE.yuv';
    outFiles.ext.pred = '_PRED.yuv';
    outFiles.mode = [writeFileFmt outFiles.ext.mode];
    outFiles.pred = [writeFileFmt outFiles.ext.pred];
end


%% YCbCr Image Conversion
LayerRec = Layer;
if(~skipYUVConv)
    LayerRec = RGB2YCbCr(LayerRec, 1);
end


%% Rescaling of the values (to keep maximum precision and avoid clipping) and corresponding adjustment of the QP.
quantizFact = 2^bitDepth-1;
QPmin = (8-bitDepth)*6;
QPmax = 51;
QPQuantLogStep = log(2)/6; %-> An increase of 6 QP steps corresponds a multiplication of the quantization step by 2 within HEVC.
minQPAdjust=-64; %Limits fixed to allow transmission of the QPAdjust value required for decoding.
maxQPAdjust=63;

minLayVal = int16(floor( min(LayerRec(:)) * quantizFact ));
QPAdjust = ceil( log(max(LayerRec(:)) - double(minLayVal)/quantizFact) / QPQuantLogStep );
if(QPAdjust < minQPAdjust)
    warning(['QP Adjustement is too low. Using minumum authorized value (QPAdjust=' num2str(minQPAdjust) ') ==> The encoding will have lower internal precision than expected.']);
    QPAdjust = minQPAdjust;
elseif(QPAdjust > maxQPAdjust)
    warning(['QP Adjustement is too high. Using maximum authorized value (QPAdjust=' num2str(minQPAdjust) ') ==> The high values in the image will be clipped.']);
    QPAdjust = maxQPAdjust;
end

QP = QPVirtual - QPAdjust; %Adjust QP with respect to the scaling applied to the layer values.

if(QP > QPmax)
    warning(['required QP (' num2str(QP) ') is higher than the maximum authorized by HEVC. Encoding using maximum QP (' num2str(QPmax) ') instead']);
    QP = QPmax;
    QPAdjust = min(QPVirtual - QP, maxQPAdjust); %Increase QPAdjust to compensate for modified QP ==> lower internal precision (not really a problem for highly compressed data).
elseif(QP < QPmin)
    warning(['required QP (' num2str(QP) ') is lower than the minimum authorized by HEVC. Encoding using minimum QP (' num2str(QPmin) ') instead']);
    QP = QPmin;
    %Here, QPAdjust isn't modified (decreased). This would clip high image values ==> Instead the image will be encoded with lower quality than requested (higher QP).
end

QPAdjust = int8(QPAdjust);


%% HEVC Encoding
[resXPad, resYPad] = writeYUV_ScaleForHEVC(outFiles.org, LayerRec, bitDepth, chromaFormat, minLayVal, QPAdjust);

nbFrames=1;
IntraPeriod = 1;
fps=25;
GenerateOptionalFiles='';
if(useRef)
    ILPOptions = [' --InputFileLDR="' refLayerFile '" --LDRLayerBitDepth=' num2str(bitDepth) ' --ILPStatus=1' ' --LambdaFactor=' num2str(lFact)];%
    if(outputModePred)
        GenerateOptionalFiles = [' --PrintMode="' outFiles.mode '" --PredFile="' outFiles.pred '"'];
    end
else
    ILPOptions = '';
end

Encod_cmd = ['"' HEVC_exe '" -c "' HEVC_cfg '" --SourceWidth=' num2str(resXPad) ' --SourceHeight=' num2str(resYPad) ...
                ' --InputFile="' outFiles.org  '" --IntraPeriod=' num2str(IntraPeriod) ' --FrameRate=' num2str(fps) ' --InputChromaFormat=' chromaFormat ...
                ' --InputBitDepth=' num2str(bitDepth) ' --InternalBitDepth=' num2str(bitDepth) ' --FramesToBeEncoded=' num2str(nbFrames) ...
                ILPOptions GenerateOptionalFiles ...
                ' --ExtendedPrecision=1  --QP=' num2str(QP) ' --ReconFile="' outFiles.rec '" --BitstreamFile="' outFiles.bin '" > "' outFiles.log '"'];

[~, errMsg] = system(Encod_cmd);
if(~exist(outFiles.bin,'file') || ~exist(outFiles.rec,'file'))
    if(exist(outFiles.log,'file')), errMsg = [errMsg 10 '--> For more detail, check log file: "' outFiles.log '"'];end
    error(['An error occured during encoding:' 10 errMsg]);
end

LayerRec = readYUV_ScaleForHEVC(outFiles.rec, resX, resY, bitDepth, chromaFormat, minLayVal, QPAdjust);

%% Inverse YCbCr Conversion
if(~skipYUVConv)
    LayerRec = YCbCr2RGB(LayerRec, 1);
end

%% Compute MSE and read bitrate
fbin = dir(outFiles.bin);
nBits = fbin.bytes*8;
if(nargout>5)
    MSE = mean((Layer(:) - LayerRec(:)).^2);
end
