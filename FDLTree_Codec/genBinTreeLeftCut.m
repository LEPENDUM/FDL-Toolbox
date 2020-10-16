%Generates a structure binary tree with leftmost node splitting:
%-A structure tree is defined so that every node's key is equal to the
%number of leaves that descend from that node (e.g. the root node's key is
%equal to the number of leaves of the tree, and each leaf node has a key
%value of 1).
%-Each parent node is split so that the left child is a leaf node (i.e.
%left child's key equal to 1 and right child's key equal to the parent's
%key minus 1).
%
%Input:
% - rootKey: Number of leaf nodes (i.e. value of the root node's key).
%
%Output:
% - tree: Generated binary tree.

function tree = genBinTreeLeftCut(rootKey)
    tree = leftSplit(BinNode(rootKey));
end

function curNode = leftSplit(curNode)
    keyCur = curNode.key;
    
    if(keyCur>1)
        keyLeft = 1;
        keyRight = curNode.key - keyLeft;
        curNode.left = BinNode(keyLeft);
        curNode.right = BinNode(keyRight);
        
        curNode.left = leftSplit(curNode.left);
        curNode.right = leftSplit(curNode.right);
    end
end