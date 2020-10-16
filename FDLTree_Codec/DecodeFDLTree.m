%Decode a FDL model encoded with Hierarchical Binary encoding method.
% Directories parameters:
% - binDir  : Directory of the encoded binary files.
% - tempDir : Temporary directory for intermediate files (these files are removed after the decoding is complete).
%             By default, tempDir is the same as binDir.
%
% Configuration parameters:
% - cfg.chromaFormat : YUV Chroma format (should be the same as the encoder's configuration).
%                      - '444' encode for full resolution chroma.
%                      - '420' chroma downsampling.
%                      - '400' only Luma.
% - cfg.bitDepth     : Encoding bitdepth in HEVC (must be the same as the encoder's configuration).
% - cfg.splitCodec   : Codec used for Split Images (should be the same as the encoder's configuration):
%                      - 'HEVC' : standard HEVC (only intra).
%                      - 'HEVC+linILP' (default) : HEVC with additional block-wise linear Inter Layer Prediction mode.
%                      - 'HEVC+nonlinILP' : HEVC with additional block-wise non-linear Inter Layer Prediction mode.
% - cfg.skipYUVConv  : (default = false) Set true to skip the YUV conversions in the decoder (must be the same as the encoder's configuration).
% - cfg.level        : If specified, decode the tree only up to this level (otherwise decode all levels).
%
% Outputs
% - FDLspat   : FDL model in the spatial domain after decoding the tree up to the maximum level (or up to cfg.level if specified).
% - Disps     : Disparity values of the FDL model associated with the layers in FDLspat.
% - nBits     : Total number of bits (including tree metadata) for decoding up to the maximum level (or up to cfg.level if specified).
% - decodTime : Decoding time in seconds (excludes the operations of reading and writing to the disk).
% - FDLTree   : Decoded FDLTree up to the maximum level (or up to cfg.level if specified).
%               See ReadFDLTree.m to retrieve the reconstructed FDL from FDLTree at a given level.
%

function [FDLspat, Disps, nBits, decodTime, FDLTree] = DecodeFDLTree(binDir, tempDir, cfg)
if(~exist(tempDir,'dir'))
    mkdir(tempDir);
end

%% Default config parameters
if(~isfield(cfg,'skipYUVConv')), cfg.skipYUVConv = false;end
if(~isfield(cfg,'splitCodec')),  cfg.splitCodec = 'HEVC+linILP';end

%% Read metada tree (treeData) with corresponding structure tree
treeMetadata_filename = fullfile(binDir,'auxData.bin');
if(~exist(treeMetadata_filename,'file'))
    error(['Can''t find tree metadata file: ' treeMetadata_filename]);
end
[resX,resY,treeData,splitTree]=CompoundLayer.readTreeMetadata(treeMetadata_filename);
numLayers = splitTree.key();

%% Decoding configuration parameters
cfg.resX=resX;
cfg.resY=resY;
cfg.binDir = binDir;
cfg.tempDir = tempDir;
cfg.fileNameFmt = 'Lvl=%d_Ids=%d-%d';
cfg.binExt = '.bin';
cfg.logExt = '.log';
cfg.decExt = '_DEC.yuv';

if(isfield(cfg,'level'))
    level = cfg.level;
else
    level = splitTree.depth;
end


%% Decode level 0 layer.
fName0 = sprintf(cfg.fileNameFmt,0,1,numLayers);
bitstreamFile0 = fullfile(cfg.binDir, [fName0 cfg.binExt]);
decLayerFile0 = fullfile(cfg.tempDir, [fName0 cfg.decExt]);
[Layer0, decodTime0] = DecodeLayer(bitstreamFile0, decLayerFile0, [], resX, resY, cfg.bitDepth, treeData.key.minLayVal, treeData.key.QPScaleAdj, cfg.chromaFormat, cfg.skipYUVConv, cfg.splitCodec);
fbin = dir(bitstreamFile0);

%% Initialize FDL tree and decode further levels recursively.
FDLTree = BinNode(CompoundLayer(0, 1, numLayers, Layer0));
FDLTree.key.decodTime = decodTime0;
FDLTree.key.nBits = fbin.bytes*8 + treeData.key.nBits; %Count the coded data and metadata of first node in the bitcount of root node.
[FDLspat,Disps,decodTime,nBits] = decodeRecur(FDLTree, treeData, splitTree, level, cfg);



%% Delete temporary folder's content
filePattern = fullfile(tempDir, ['*' cfg.decExt]);
f = dir(filePattern);
for i=1:length(f)
  fileName = fullfile(tempDir, f(i).name);
  delete(fileName);
end
logFile = fullfile(tempDir, 'dec.log');
if(exist(logFile,'file'))
    delete(logFile);
end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Recursive tree decoding function
function [FDLspat,Disps,decodTime,nBits] = decodeRecur(curNodeFDL, curNodeAux, curNodeSplit, tgtLvl, cfg)
    if(~isempty(curNodeAux))
        curLvl = curNodeFDL.key.level;
        curNodeFDL.key.disp = curNodeAux.key.disp;
        if(tgtLvl==curLvl || curNodeAux.isLeaf())
            FDLspat = curNodeFDL.key.Layer;
            Disps = curNodeAux.key.disp;
            decodTime = curNodeFDL.key.decodTime;
            nBits = curNodeFDL.key.nBits;
        else
            
            minLayVal  = curNodeAux.left.key.minLayVal;
            QPScaleAdj = curNodeAux.left.key.QPScaleAdj;
            
            numL = curNodeSplit.left.key;
            numR = curNodeSplit.right.key;
            idRangeL = [curNodeFDL.key.idStart, curNodeFDL.key.idStart+numL-1];
            idRangeR = [curNodeFDL.key.idStart+numL, curNodeFDL.key.idStart+numL+numR-1];
            idRange = [curNodeFDL.key.idStart, curNodeFDL.key.idEnd];
            
            %Set filename variables
            fNameL = sprintf(cfg.fileNameFmt,curLvl+1,idRangeL(1),idRangeL(2));
            fNameR = sprintf(cfg.fileNameFmt,curLvl+1,idRangeR(1),idRangeR(2));
            decLayerFileL = fullfile(cfg.tempDir, [fNameL cfg.decExt]);
            decLayerFileR = fullfile(cfg.tempDir, [fNameR cfg.decExt]);
            refLayerFile = fullfile(cfg.tempDir, [sprintf(cfg.fileNameFmt,curLvl,idRange(1),idRange(2)) cfg.decExt]);
            bitstreamFile = fullfile(cfg.binDir, [fNameL cfg.binExt]);
            decLayerFile = decLayerFileL;
            
            %Create the children nodes.
            curNodeFDL.left = BinNode(CompoundLayer(curLvl+1, idRangeL(1), idRangeL(2)));
            curNodeFDL.right = BinNode(CompoundLayer(curLvl+1, idRangeR(1), idRangeR(2)));
            
            
            %Decode the Layers:
            % -> Decode image data.
            [LayerDiffDec, decodTimeLayer] = DecodeLayer(bitstreamFile, decLayerFile, refLayerFile, cfg.resX, cfg.resY, cfg.bitDepth, minLayVal, QPScaleAdj, cfg.chromaFormat, cfg.skipYUVConv, cfg.splitCodec);
            
            % -> Reconstruct the two sub-layers from the parent layer and the decoded image.
            start = tic;
            ratio_rl = numL/numR;%numR / numL;%
            curNodeFDL.left.key.Layer = curNodeFDL.key.Layer  * ratio_rl/(1+ratio_rl) + LayerDiffDec;
            curNodeFDL.right.key.Layer = curNodeFDL.key.Layer * 1/(1+ratio_rl)        - LayerDiffDec;
            curNodeFDL.right.key.decodTime = toc(start);    %right node stores final reconstruction time.
            curNodeFDL.left.key.decodTime = decodTimeLayer; %left node stores HEVC decoding time.

            %Read the encoding log to retrieve number of bits data (excluding HEVC header).
            fbin = dir(bitstreamFile);
            curNodeFDL.left.key.nBits = fbin.bytes*8 + curNodeAux.left.key.nBits; %Count the coded data and metadata in the bitcount of left node.
            curNodeFDL.right.key.nBits = 0;

            %Write file of the reconstructed compound layers for later prediction of their children layers.
            if(~(curNodeAux.left.isLeaf))
                LayerWrite = curNodeFDL.left.key.Layer;
                LayerWrite = (LayerWrite - min(LayerWrite(:))) / (max(LayerWrite(:)) - min(LayerWrite(:)));
                if(~cfg.skipYUVConv), LayerWrite = RGB2YCbCr(LayerWrite,1); end
                writeYUV_ScaleForHEVC(decLayerFileL, LayerWrite, cfg.bitDepth, cfg.chromaFormat, 0, 0);
            end
            if(~(curNodeAux.right.isLeaf))
                LayerWrite = curNodeFDL.right.key.Layer;
                LayerWrite = (LayerWrite - min(LayerWrite(:))) / (max(LayerWrite(:)) - min(LayerWrite(:)));
                if(~cfg.skipYUVConv), LayerWrite = RGB2YCbCr(LayerWrite,1); end
                writeYUV_ScaleForHEVC(decLayerFileR, LayerWrite, cfg.bitDepth, cfg.chromaFormat, 0, 0);
            end
            
            %Decode next levels recursively
            if(~isempty(curNodeAux.left))
                [FDLspatL,DispsL,decodTimeL,nBitsL] = decodeRecur(curNodeFDL.left, curNodeAux.left, curNodeSplit.left, tgtLvl, cfg);
            end
            if(~isempty(curNodeAux.right))
                [FDLspatR,DispsR,decodTimeR,nBitsR] = decodeRecur(curNodeFDL.right, curNodeAux.right, curNodeSplit.right, tgtLvl, cfg);
            end
            
            %Concatenate result layers and disparity values from target level.
            FDLspat = cat(4,FDLspatL,FDLspatR);
            Disps = [DispsL, DispsR];
            decodTime = decodTimeL + decodTimeR + curNodeFDL.key.decodTime;
            nBits = nBitsL + nBitsR + curNodeFDL.key.nBits;
            
        end
    else
        FDLspat=[];
        Disps=[];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Layer decoding function
function [LayerDec, decodTime] = DecodeLayer(bitstreamFile, decLayerFile, refLayerFile, resX, resY, bitDepth, minLayVal, QPAdjust, chromaFormat, skipYUVConv, splitCodec)

if(isempty(refLayerFile))
    HEVC_exe = [HEVC_ILP_path 'TAppDecoder_HM-12.1+RExt-4.2.exe'];
    refCmd = '';
else
    refCmd = [' --InputFileLDR="' refLayerFile '" --LDRLayerBitDepth=' num2str(bitDepth)];
    switch splitCodec
        case 'HEVC'
            HEVC_exe = [HEVC_ILP_path 'TAppDecoder_HM-12.1+RExt-4.2.exe'];
            refCmd = '';
        case 'HEVC+linILP'
            HEVC_exe = [HEVC_ILP_path 'TAppDecoder_ILP_lin_tr.exe'];
        case 'HEVC+nonlinILP'
            HEVC_exe = [HEVC_ILP_path 'TAppDecoder_ILP_3lin_Ext_Hybrid.exe'];
        otherwise
            error(['Unknown decoding type ''' splitCodec ''' specified in the splitCodec configuration field.']);
    end
end

logFile = fullfile(fileparts(decLayerFile), 'dec.log');

Decod_cmd = ['"' HEVC_exe '" -b "' bitstreamFile '"' refCmd ' --ReconFile="' decLayerFile '" > "' logFile '"'];
[~, errMsg] = system(Decod_cmd);
if(~exist(decLayerFile,'file'))
    if(exist(logFile,'file')), errMsg = [errMsg 10 '--> For more detail, check log file: "' logFile '"'];end
    error(['An error occured during decoding:' 10 errMsg]);
end

LayerDec = readYUV_ScaleForHEVC(decLayerFile, resX, resY, bitDepth, chromaFormat, minLayVal, QPAdjust);
decodTime = parse_HEVC_dec_log(logFile,1);

%% Inverse YCbCr Conversion
if(~skipYUVConv)
    LayerDec = YCbCr2RGB(LayerDec, 1);
end

end
