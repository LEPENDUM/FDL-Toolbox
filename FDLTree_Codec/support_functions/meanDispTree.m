%Generates a binary tree with specified structure and with node's keys
%corresponding to disparity values. The disparity at each node is computed
%as the mean of the disparities of the leaf nodes that descend from it.
%
%Input:
% - Disps: list of disparity values of leaf nodes (from left to right).
% - StructureTree: Structure binary tree (e.g. see genBinTreeMidCut or
% genBinTreeLeftCut). The generated tree will have the same structure.
%
%Output:
% - DispTree: Generated binary tree.

function DispTree = meanDispTree(Disps,StructureTree)
    DispTree = meanDispTreeRecur(BinNode([]),Disps,StructureTree);
end

function curNode = meanDispTreeRecur(curNode,curDisps,curStructTree)

    curNode.key = mean(curDisps);
    if(~isempty(curStructTree.left))
        split = curStructTree.left.key;
        curNode.left = meanDispTreeRecur(BinNode([]), curDisps(1:split), curStructTree.left);
    end

    if(~isempty(curStructTree.right))
        numElts = curStructTree.right.key;
        if(isempty(curStructTree.left))
            split = 0;
        else
            split = curStructTree.left.key;
        end
        curNode.right = meanDispTreeRecur(BinNode([]), curDisps(split+1:split+numElts), curStructTree.right);
    end

end