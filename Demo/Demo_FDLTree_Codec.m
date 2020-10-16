%-------------------------------------------------------------------------%
%                   Demo script for the FDL Tree streaming
%
%   See the paper: M. Le Pendu, C. Ozcinar and A. Smolic, "Hierarchical
%   Fourier Disparity Layer Transmission for Light Field Streaming", ICIP
%   2020.
%
% -> This script generates a FDL, encode it as a FDLTree, and decode the
% tree. The FDLs obtained by decoding the tree up to different levels are
% then displayed in the FDL rendering application.
%
% -> The difference in encoding quality between different levels of
% decoding can be seen by moving to external views or by refocusing using
% 'aperture radius' and 'refocus' sliders (the quality of the all-in focus
% central view does not depend on the number of levels decoded).
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construction of a FDL (uncompressed). %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Demo_FDL_Simple;
figInfo = gcf; figNum = figInfo.Number;
disp([10 '---------------------------' 10 ...
         'Figure ' num2str(figNum) ': Uncompressed FDL' 10 ...
         '---------------------------' 10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Encoding configuration  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QP = 22;

cfg = struct(...
    'QPBase', QP,...
    'bitDepth', 10,...
    'chromaFormat','420',...
    'splitCodec','HEVC+linILP');
resDir = fullfile(LFDir,[LFName '_QP' num2str(QP)]);
resBinDir = fullfile(resDir,'bin');
tempDir = fullfile(LFDir,'temp');


%% Encode the FDL using FDL Tree codec.
FDLspat = ifftImgs(FDL);
FDLTree_rec = EncodeFDLTree(FDLspat, D, genBinTreeMidCut(numLayers), cfg, resDir);


%% Decode the FDL Tree from bitstream.
[FDLspat_decFull, D_decFull, nBits_Full, decodTime_Full, FDLTree_decFull] = DecodeFDLTree(resBinDir,tempDir, cfg);

% Start Rendering Application -> fully decoded
RenderAppMain(fftImgs(FDLspat_decFull), [fullResY fullResX], bp, D_decFull,U,V, [],{linearizeInput,gammaOffset}, useGPU);

figInfo = gcf; figNum = figInfo.Number;
disp([10 '---------------------------------' 10 ...
'Figure ' num2str(figNum) ': Decoded FDLTree (Full)' 10 ...
' number of bits: ' num2str(nBits_Full) 10 ...
' number of layers: ' num2str(length(D_decFull)) 10 ...
'----------------------------------' 10]);


%% Read decoded data only up to a given level of the tree (same result could be obtained by decoding directly with DecodeFDLTree using the parameter cfg.level).
lvl_max = 2;
[FDLspat_decLvl, D_decLvl, nBits_decLvl, decodTime_decLvl] = ReadFDLTree(FDLTree_decFull,lvl_max);

% Start Rendering Application -> decoded up to level 'lvl_max'
RenderAppMain(fftImgs(FDLspat_decLvl), [fullResY fullResX], bp, D_decLvl,U,V, [],{linearizeInput,gammaOffset}, useGPU);

figInfo = gcf; figNum = figInfo.Number;
disp([10 '------------------------------------' 10 ...
'Figure ' num2str(figNum) ': Decoded FDLTree (level ' num2str(lvl_max) ')' 10 ...
' number of bits: ' num2str(nBits_decLvl) 10 ...
' number of layers: ' num2str(length(D_decLvl)) 10 ...
'------------------------------------' 10]);


%% Read decoded data only up to a given level of the tree (same result could be obtained by decoding directly with DecodeFDLTree using the parameter cfg.level).
lvl_max = 1;
[FDLspat_decLvl, D_decLvl, nBits_decLvl, decodTime_decLvl] = ReadFDLTree(FDLTree_decFull,lvl_max);

% Start Rendering Application -> decoded up to level 'lvl_max'
RenderAppMain(fftImgs(FDLspat_decLvl), [fullResY fullResX], bp, D_decLvl,U,V, [],{linearizeInput,gammaOffset}, useGPU);

figInfo = gcf; figNum = figInfo.Number;
disp([10 '------------------------------------' 10 ...
'Figure ' num2str(figNum) ': Decoded FDLTree (level ' num2str(lvl_max) ')' 10 ...
' number of bits: ' num2str(nBits_decLvl) 10 ...
' number of layers: ' num2str(length(D_decLvl)) 10 ...
'------------------------------------' 10]);
