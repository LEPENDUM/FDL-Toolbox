%Class for storing a Compound Layer and associated metadata.
%A FDLTree is a tree of CompoundLayers (i.e. BinNode object with key value
%of type CompoundLayer)
%
% This class contains the static functions CompoundLayer.writeTreeMetadata
% and CompoundLayer.readTreeMetadata for writting and reading the metadata
% of a FDLTree required for decoding:
%   - Horizontal and vertical resolution of the image : 2*16-bit unsigned integers.
%   - Structure of the FDL tree: number of bits = 2*number_of_leaves-1.
%   - Required metadata for the root node and each left node: number of bits = number_of_leaves*(32+16+7).
%       - disp: single (32bit).
%       - QPScaleAdj: 7-bit signed integer.
%       - minLayVal: 16-bit signed integer.

classdef CompoundLayer
    
    properties
        level=0
        idStart=1
        idEnd=1
        Layer=[]
        disp=[]
        nBits=0
        MSE=0
        minLayVal=0
        QPScaleAdj=0
        decodTime=0
    end
    
    methods
        
        function obj = CompoundLayer(level, idStart, idEnd, Layer, disp, nBits, MSE, minLayVal, QPScaleAdj)
            if(nargin>0)
                obj.level=level;
            end
            if(nargin>1)
                obj.idStart = idStart;
            end
            if(nargin>2)
                obj.idEnd = idEnd;
            end
            if(nargin>3)
                obj.Layer = Layer;
            end
            if(nargin>4)
                obj.disp = disp;
            end
            if(nargin>5)
                obj.nBits = nBits;
            end
            if(nargin>6)
                obj.MSE = MSE;
            end
            if(nargin>7)
                obj.minLayVal = minLayVal;
            end
            if(nargin>8)
                obj.QPScaleAdj = QPScaleAdj;
            end
        end
        
    end
    
    methods( Static )
        %Write metadata of a FDLTree.
        function writeTreeMetadata(FDLTree, filename)
            
            resX = size(FDLTree.key.Layer,2);
            resY = size(FDLTree.key.Layer,1);
            if(~isequal(resX, uint16(resX))), error('Image horizontal resolution must fit in uint16 format');end
            if(~isequal(resX, uint16(resX))), error('Image vertical resolution must fit in uint16 format');end
            
            bitstream = BitStream();
            %add header metadata to bitstream
            bitstream.add_unsignedN(uint16(resX),16);
            bitstream.add_unsignedN(uint16(resY),16);
            FDLTree.structure2bitstream(bitstream); %add structure of the tree to the bitstream.
            FDLTree.key.nBits = FDLTree.key.nBits + bitstream.addedCount();%Count header bits in the root node.
            %add tree metadata to bitstream
            CompoundLayer.treeMetadata2bitstream(FDLTree,bitstream); %add metadata of the tree nodes to the bitstream.
            
            fid = fopen(filename,'w');
            try
                bitstream.writeToFile(fid); 
            %close file even if an error is raised during execution.
                fclose(fid);
            catch exception
                fclose(fid);
                rethrow(exception);
            end
        end
        
        %Read metadata of a FDLTree from file.
        function [resX,resY,FDLTree_meta,structureTree]=readTreeMetadata(filename)
            bitstream = BitStream();
            fid = fopen(filename);
            try
                bitstream.addFromFile(fid);
            %close file even if an error is raised during execution.
                fclose(fid);
            catch exception
                fclose(fid);
                rethrow(exception);
            end
            
            %Read header metadata
            resX = bitstream.read_unsignedN(16);
            resY = bitstream.read_unsignedN(16);
            structureTree = BinNode.bitstream2tree(bitstream); %Read structure of the tree from the bitstream.
            structureTree = structureTree.genStructureTree();  %Convert to structure tree (i.e. keys equal to number of descending leaves).
            nBitsHeader = bitstream.readCount();
            
            %Read tree metadata
            FDLTree_meta = CompoundLayer.bitstream2treeMetadata(bitstream, structureTree); %Read tree metadata.
            FDLTree_meta.key.nBits = FDLTree_meta.key.nBits + nBitsHeader;%Count header bits in the root node.
        end
        
    end
    
    methods( Access = private, Static )
        %Add metadata of a single node to the bitstream.
        function metadata2bitstream(metadata,bitstream)
            bitstream.add_single(metadata.disp);                 %disp -> single: 
            bitstream.add_signedN(int8(metadata.QPScaleAdj),7);  %QPScaleAdj -> signed 7 bit integer 
            bitstream.add_signedN(metadata.minLayVal,16);        %minLayVal -> int16
        end
        %Read metadata of a single node from the bitstream.
        function metadata = bitstream2metadata(bitstream)
            metadata = struct(...
                'disp',       bitstream.read_single(),...
                'QPScaleAdj', bitstream.read_signedN(7),...
                'minLayVal',  bitstream.read_signedN(16));
        end
        
        %Binarize metadata of the tree nodes in breadth-first order.
        function treeMetadata2bitstream(FDLTree, bitstream)
            if(~isempty(FDLTree))
                bitstream.resetAddedCount();
                CompoundLayer.metadata2bitstream(FDLTree.key, bitstream);
                FDLTree.key.nBits = FDLTree.key.nBits + bitstream.addedCount();
            end
            queue = FDLTree;
            while ~isempty(queue)
                curNode = queue(1);
                queue(1)=[];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(~curNode.isLeaf())
                    bitstream.resetAddedCount();
                    CompoundLayer.metadata2bitstream(curNode.left.key, bitstream);
                    curNode.left.key.nBits = curNode.left.key.nBits + bitstream.addedCount();
                    
                    %Compute right node disparity from parent and left node to avoid transmitting additional metadata.
                    %If right node disparity value already exists in the tree and is different from the computed value, an error is raised.
                    nl = 1 + curNode.left.key.idEnd - curNode.left.key.idStart;
                    nr = 1 + curNode.right.key.idEnd - curNode.right.key.idStart;
                    rightDisp_rec = single(((nr+nl)/nr)*(double(curNode.key.disp)) - (nl/nr)*double(curNode.left.key.disp));
                    if(isempty(curNode.right.key.disp))
                        curNode.right.key.disp = rightDisp_rec;
                    else
                        if(rightDisp_rec ~= curNode.right.key.disp), error('Right node disparity in FDL Tree does not match reconstructed disparity value');end
                    end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    queue(end+1) = curNode.left;
                    queue(end+1) = curNode.right;
                end
            end
        end
        
        %Read Tree metadata from binary data, given the structure tree as additional input.
        function FDLTree = bitstream2treeMetadata(bitstream, structTree)
            if(isempty(structTree) || ~isequal(class(structTree),'BinNode'))
                FDLTree = []; return;
            else
                FDLTree = BinNode(CompoundLayer());
                
                bitstream.resetReadCount();
                node_meta = CompoundLayer.bitstream2metadata(bitstream);
                
                FDLTree.key.disp = node_meta.disp;
                FDLTree.key.QPScaleAdj = node_meta.QPScaleAdj;
                FDLTree.key.minLayVal = node_meta.minLayVal;
                FDLTree.key.nBits = bitstream.readCount();
            end
            
            queueStruct = structTree;
            queueFDL = FDLTree;
            while ~bitstream.isempty() && ~isempty(queueStruct)
                curNodeStruct = queueStruct(1);
                queueStruct(1)=[];
                curNode = queueFDL(1);
                queueFDL(1)=[];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(~curNodeStruct.isLeaf())
                    nl = curNodeStruct.left.key;
                    nr = curNodeStruct.right.key;
                    
                    bitstream.resetReadCount();
                    node_meta = CompoundLayer.bitstream2metadata(bitstream);
                    
                    curNode.left = BinNode(CompoundLayer());
                    curNode.left.key.disp = node_meta.disp;
                    curNode.left.key.QPScaleAdj = node_meta.QPScaleAdj;
                    curNode.left.key.minLayVal = node_meta.minLayVal;
                    curNode.left.key.nBits = bitstream.readCount();
                    
                    curNode.right = BinNode(CompoundLayer());
                    curNode.right.key.disp = single(((nr+nl)/nr)*(double(curNode.key.disp)) - (nl/nr)*double(node_meta.disp)); %find right node disparity from parent and left nodes.
                    curNode.right.key.minLayVal = int16(0);
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    queueStruct(end+1) = curNodeStruct.left;
                    queueStruct(end+1) = curNodeStruct.right;
                    queueFDL(end+1) = curNode.left;
                    queueFDL(end+1) = curNode.right;
                end
            end
            if(~isempty(queueStruct))
                error('Bitstream incompatible with metadata tree structure (may be a truncated bitstream)');
            end
        end
    end
    
end