%Read a given level of a FDL tree to convert it into an array structure.
%This function assumes that each node contains a reconstructed FDL layer
%(i.e. NOT the intermediate split image).
%
%Inputs:
% - FDLTree: FDL stored in a tree structure. Each node should contain a
% reconstructed FDL layer (i.e. NOT the intermediate split image).
% - level : level of the tree to read.
%
% Outputs
% - FDLspat   : FDL array (spatial domain).
% - Disps     : Disparity values of the FDL model.
% - nBits     : Number of bits up to the read level (assumes number of bits
% are given for each node).
% - decodTime : Sum of decoding times of each nodes up to the read level
% (assumes decoding times are given for each node).

function [FDLspat, Disps, nBits, decodTime] = ReadFDLTree(FDLTree,level)

numLayers = FDLTree.numLeavesUpToLevel(level);

FDLSpat0_rec = FDLTree.key.Layer;

imgSize = [size(FDLSpat0_rec,1) size(FDLSpat0_rec,2)];
nChan = size(FDLSpat0_rec, 3);

%Initialize Data Structures
Disps = zeros(1, numLayers);
FDLspat = zeros([imgSize nChan numLayers]);

[~,FDLspat,Disps,nBits,decodTime] = readLevelRecur(FDLTree,level,0,1,FDLspat,Disps);

end

function [curId,FDLspat,Disps,nBits,decodTime] = readLevelRecur(curNode,tgtLvl,curLvl,curId,FDLspat,Disps)
    if(~isempty(curNode))        
        if(tgtLvl==curLvl || curNode.isLeaf())
            FDLspat(:,:,:,curId) = curNode.key.Layer;
            Disps(curId) = curNode.key.disp;
            nBits = curNode.key.nBits;
            decodTime = curNode.key.decodTime;
            curId = curId+1;
        else
            if(~isempty(curNode.left))
                [curId,FDLspat,Disps,nBitsL,decodTimeL] = readLevelRecur(curNode.left, tgtLvl, curLvl+1, curId, FDLspat, Disps);
            end
            if(~isempty(curNode.right))
                [curId,FDLspat,Disps,nBitsR,decodTimeR] = readLevelRecur(curNode.right, tgtLvl, curLvl+1, curId, FDLspat, Disps);
            end
            nBits = nBitsL + nBitsR + curNode.key.nBits;
            decodTime = decodTimeL + decodTimeR + curNode.key.decodTime;
        end
    else
        nBits=0;
        decodTime=0;
    end
end