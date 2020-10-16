%Encode a FDL model with Hierarchical Binary encoding method.
% Input Data
% - FDLspat: FDL model converted in the spatial domain (dimensions : resY, resX, #Channels, #Layers).
% - Disps: list of disparity values of the FDL (one value per layer).
%
% Encoding Configuration:
% - splitTree: Structure Tree of the compound layer hierarchy (can be generated with e.g. genBinTreeMidCut, genBinTreeLeftCut).
% - cfg.QPBase : Base QP of the HEVC encoding.
% - cfg.chromaFormat : YUV Chroma format:
%                      - '444' encode for full resolution chroma.
%                      - '420' chroma downsampling.
%                      - '400' only Luma.
% - cfg.bitDepth     : Encoding bitdepth in HEVC.
% - cfg.splitCodec   : Codec used for Split Images:
%                      - 'HEVC' : standard HEVC (only intra).
%                      - 'HEVC+linILP' (default) : HEVC with additional block-wise linear Inter Layer Prediction mode.
%                      - 'HEVC+nonlinILP' : HEVC with additional block-wise non-linear Inter Layer Prediction mode.
% - cfg.skipYUVConv  : (default = false) Set true to skip the YUV conversions in the encoder (e.g. if the input FDL is already given in YUV colorspace).
% - cfg.lFact        : (default=1) Factor to multiply the lagrangian parameter lambda used in HEVC for rate distortion optimisation (only used if cfg.splitCodec is 'HEVC+linILP' or 'HEVC+nonlinILP').
% - cfg.dQPLvl       : (default = 0) deltaQP per Level => Value to add to the base QP at each level of the encoding (use negative values to have lower QPs for higher levels).
% - cfg.outputModePred : (default=false) Set true to generate prediction images and mode decision images (only used if cfg.splitCodec is 'HEVC+linILP' or 'HEVC+nonlinILP').
%
% Output Directories:
% - resDir  : Result directory to save bitstream and log files.
% - tempDir : Temporary directory for intermediate files (these files are removed after the encoding).
%             By default tempDir in the same as resDir.
%
% Output
% - FDLTree : Tree structure containing the reconstructed layers at each node.
%             See ReadFDLTree.m to retrieve the reconstructed FDL from FDLTree at a given level.
% The function generates subdirectories 'bin' and 'log' in the resDir
% directory containing respectively the binary files of the full encoded
% tree and the corresponding HEVC log files.

function FDLTree = EncodeFDLTree(FDLspat, Disps, splitTree, cfg, resDir, tempDir)
numLayers = size(FDLspat,4);
if(length(Disps)~=numLayers)
    error('The number of elements in the list of disparity values should be equal to the number of layers.');
end

%% Default config parameters
if(~isfield(cfg,'dQPLvl')), cfg.dQPLvl = 0;end
if(~isfield(cfg,'skipYUVConv')), cfg.skipYUVConv = false;end
if(~isfield(cfg,'splitCodec')),  cfg.splitCodec = 'HEVC+linILP';end


%% Output Directory and filename variables
resDir = fullfile(resDir,'/');
if(~exist('tempDir','var') || isempty(tempDir))
    tempDir = resDir;
else
    tempDir = fullfile(tempDir,'/');
end
bitstreamDir = fullfile(resDir, 'bin/');
logDir       = fullfile(resDir, 'log/');

tmpFileFmt  = strrep([tempDir 'Lvl=%d_Ids=%d-%d_QPb=' num2str(cfg.QPBase)],'\','\\');
resFileNamesFmt = 'Lvl=%d_Ids=%d-%d';

%Create all directories
if(~exist(resDir,'dir'))
    mkdir(resDir);
end
if(~exist(tempDir,'dir'))
    mkdir(tempDir);
end
if(~exist(logDir,'dir'))
    mkdir(logDir);
end
if(~exist(bitstreamDir,'dir'))
    mkdir(bitstreamDir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Encode initial layer (sum of all layers = rendered central view).

FDLspat_org = sum(FDLspat,4);
disp0 = single(mean(double(Disps)));

%Encoding (without reference layer).
fprintf('Encoding Level 0...');
recFileName0 = sprintf(tmpFileFmt,0,1,numLayers);
[FDLSpat0_rec, nBits0, fNames0, min0, QPAdj0, MSE0] = EncodeLayer(FDLspat_org, recFileName0, [], cfg.QPBase, cfg);
%Move generated files from temporary directory to final results directory
movefile(fNames0.log, [logDir sprintf(resFileNamesFmt,0,1,numLayers) fNames0.ext.log]);
movefile(fNames0.bin, [bitstreamDir sprintf(resFileNamesFmt,0,1,numLayers) fNames0.ext.bin]);
delete(fNames0.org);
%Initialize FDL tree
FDLTree = BinNode(CompoundLayer(0, 1, numLayers, FDLSpat0_rec, disp0, nBits0, MSE0, min0, QPAdj0));
fprintf('\b\b\b -> Done.\n');

%% Encode FDL hierarchically with successive split.
prevSplitLevel=-inf;
if(splitTree.isLeaf)
    SplitQueue=[];
else
    SplitQueue = splitTree;
end
FDLQueue = FDLTree;
while ~isempty(SplitQueue)
    curSplitNode = SplitQueue(1);
    SplitQueue(1)=[];
    
    curFDLNode = FDLQueue(1);
    FDLQueue(1)=[];
    
    %case with at most one child => no need to split => nothing to encode.
    if(curSplitNode.numChildren<2)
        if(isempty(SplitQueue))
           fprintf('\b\b\b -> Done.\n');
        end
        continue;
    end
    
    curSplitLevel = curFDLNode.key.level + 1;
    if(curSplitLevel > prevSplitLevel)
        if(prevSplitLevel>=0)
            fprintf('\b\b\b -> Done.\n');
        end
        countInLevel = 1;
        fprintf('Encoding Level %d - SubLayer: ', curSplitLevel);
    else
        countInLevel = countInLevel+1;
    end
    if(countInLevel>1)
        fprintf('\b\b\b / ');
    end
    fprintf('%d...',countInLevel);
    prevSplitLevel = curSplitLevel;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numL = curSplitNode.left.key;
    numR = curSplitNode.right.key;
    NodeL = CompoundLayer(curFDLNode.key.level+1, curFDLNode.key.idStart,      curFDLNode.key.idStart+numL-1);
    NodeR = CompoundLayer(curFDLNode.key.level+1, curFDLNode.key.idStart+numL, curFDLNode.key.idStart+numL+numR-1);
    refFileName = sprintf(tmpFileFmt, curFDLNode.key.level, curFDLNode.key.idStart, curFDLNode.key.idEnd);
    encFileName = sprintf(tmpFileFmt, NodeL.level, curFDLNode.key.idStart, curFDLNode.key.idEnd);
    recFileNameL = sprintf(tmpFileFmt, NodeL.level, NodeL.idStart, NodeL.idEnd);
    recFileNameR = sprintf(tmpFileFmt, NodeR.level, NodeR.idStart, NodeR.idEnd);
    resFileName = sprintf(resFileNamesFmt, NodeL.level, NodeL.idStart, NodeL.idEnd);
    
    %Determine QP parameter for the level.
    QPLevel = cfg.QPBase + cfg.dQPLvl*curSplitLevel;
    
    %Encoding of the Sub-layers...
    % -> Encode the difference image.
    FDLspatLeft_org = sum(FDLspat(:,:,:,NodeL.idStart:NodeL.idEnd), 4);
    FDLspatRight_org = sum(FDLspat(:,:,:,NodeR.idStart:NodeR.idEnd), 4);
    ratio_rl = numL/numR;
    FDLspatDiff_org = ( FDLspatLeft_org - ratio_rl * FDLspatRight_org ) / (1+ratio_rl);%difference image containing node split information.
    [FDLspatDiff_rec, nBitsDiff, fNames, minDiff, QPAdjDiff] = EncodeLayer(FDLspatDiff_org, encFileName, refFileName, QPLevel, cfg);
    
    % -> Reconstruct the two sub-layers from the parent layer and the difference image.
    NodeL.Layer = curFDLNode.key.Layer * ratio_rl/(1+ratio_rl) + FDLspatDiff_rec;
    NodeR.Layer = curFDLNode.key.Layer * 1/(1+ratio_rl)        - FDLspatDiff_rec;
    NodeL.MSE = mean((NodeL.Layer(:) - FDLspatLeft_org(:)).^2);
    NodeR.MSE = mean((NodeR.Layer(:) - FDLspatRight_org(:)).^2);
    
    % -> Set nodes auxiliary data : only needs to be stored for 1 of the children node (left child chosen here arbitrarily).
    NodeL.nBits = nBitsDiff;
    NodeL.minLayVal = minDiff;
    NodeL.QPScaleAdj = QPAdjDiff;
    NodeL.disp = single(mean(double(Disps(NodeL.idStart : NodeL.idEnd))));
    %Right node disparity is not set here because it is not transmitted: Instead it will be computed from parent and left node when writting/reading the tree metadata.
    %Note: this may cause negligible loss in the transmission of disparity values due to floating point precision errors.
    
    %Write reconstructed layer yuv files to perform prediction of lower levels.
    if(~(curSplitNode.left.isLeaf))
        LayerWrite = NodeL.Layer;
        LayerWrite = (LayerWrite - min(LayerWrite(:))) / (max(LayerWrite(:)) - min(LayerWrite(:)));
        if(~cfg.skipYUVConv), LayerWrite = RGB2YCbCr(LayerWrite,1); end
        writeYUV_ScaleForHEVC([recFileNameL fNames.ext.rec], LayerWrite, cfg.bitDepth, cfg.chromaFormat, 0, 0);
    end
    if(~(curSplitNode.right.isLeaf))
        LayerWrite = NodeR.Layer;
        LayerWrite = (LayerWrite - min(LayerWrite(:))) / (max(LayerWrite(:)) - min(LayerWrite(:)));
        if(~cfg.skipYUVConv), LayerWrite = RGB2YCbCr(LayerWrite,1); end
        writeYUV_ScaleForHEVC([recFileNameR fNames.ext.rec], LayerWrite, cfg.bitDepth, cfg.chromaFormat, 0, 0);
    end
    
    %Save log and bitream files to the result directory and remove temporary yuv files.
    movefile(fNames.log, [logDir resFileName fNames.ext.log]);
    movefile(fNames.bin, [bitstreamDir resFileName fNames.ext.bin]);
    delete(fNames.rec);
    delete(fNames.org);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curFDLNode.left = BinNode(NodeL);
    curFDLNode.right = BinNode(NodeR);
    
    SplitQueue(end+1) = curSplitNode.left;
    SplitQueue(end+1) = curSplitNode.right;
    FDLQueue(end+1) = curFDLNode.left;
    FDLQueue(end+1) = curFDLNode.right;
    
end


%% save tree metadata file.
CompoundLayer.writeTreeMetadata(FDLTree,[bitstreamDir 'auxData.bin']);

%% Delete temporary folder's content
filePattern = fullfile(tempDir, ['*' fNames0.ext.rec]);
f = dir(filePattern);
for i=1:length(f)
  fileName = fullfile(tempDir, f(i).name);
  delete(fileName);
end

end