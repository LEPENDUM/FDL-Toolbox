%Generates a structure binary tree with node splitting in the middle:
%-A structure tree is defined so that every node's key is equal to the
%number of leaves that descend from that node (e.g. the root node's key is
%equal to the number of leaves of the tree, and each leaf node has a key
%value of 1).
%-Each parent node is split so that both children have equal keys (or with
%a difference of 1 if the parent's key is not even).
%-If a parent's key is not even, the left child's key is taken lower than
%the right child's key by 1.
%
%Input:
% - rootKey: Number of leaf nodes (i.e. value of the root node's key).
%
%Output:
% - tree: Generated binary tree.

function tree = genBinTreeMidCut(rootKey)
    tree = midSplit(BinNode(rootKey));
end

function curNode = midSplit(curNode)
    keyCur = curNode.key;
    
    if(keyCur>1)
        keyLeft = floor(curNode.key/2);
        keyRight = curNode.key - keyLeft;
        curNode.left = BinNode(keyLeft);
        curNode.right = BinNode(keyRight);
        
        curNode.left = midSplit(curNode.left);
        curNode.right = midSplit(curNode.right);
    end
end