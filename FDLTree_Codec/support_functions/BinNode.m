%Class defining the data structure of binary trees.
%A BinNode object correponds to a node of a binary tree containing a key
%and two children BinNodes (left and right). A leaf node has both it's
%children empty.
%--------------------------------------------------------------------------
%
% The BinNode object constructor has the following arguments:
% -key: key of the node (must be specified).
% -left: (optional) BinNode object corresponding to the left child (set to
% empty if not specified).
% -right: (optional) BinNode object corresponding to the right child (set
% to empty if not specified).
%
%--------------------------------------------------------------------------
%
% The following public methods are  available in the BinNode class. They
% assume the caller object 'obj' is the root of a tree:
%
% - obj.isLeaf():
%       Returns true if the node is a leaf, false otherwise.
% 
% - obj.numChildren():
%       Returns the number of children nodes (0, 1 or 2).
%
% - obj.depth():
%       Returns the depth of the tree (start count from zero at root node).
%
% - obj.numLeaves():
%       Returns the number of leaves of the tree.
%
% - obj.numNodesAtLevel(level):
%       Returns the number of nodes in a given level (level number starts
%       from zero at the root node).
%
% - obj.treeUpToLevel(level):
%       Returns the the tree obtained by retaining the nodes only up to the
%       specified level (level number starts from zero at the root node).
%       All the nodes of higher level are removed from the returned tree.
%       If 'level' is higher than the tree depth, the tree depth is used.
%
% - obj.numLeavesUpToLevel(level):
%       Returns the number of leaves of the tree obtained by retaining the
%       nodes only up to the specified level (level number starts from zero
%       at the root node).
%       Always equal to obj.treeUpToLevel(level).numLeaves()
%       If 'level' is higher than the tree depth, the tree depth is used.
%
% - obj.listLeavesUpToLevel(level):
%       Returns a cell array of keys corresponding to the leaves of the
%       tree obtained by retaining the nodes only up to the specified level
%       (level number starts from zero at the root node).
%       If 'level' is higher than the tree depth, the tree depth is used.
%
% - obj.structureDiff(otherTree):
%       Compare the structure of the tree obj with the tree 'otherTree'.
%       Returns an integer flag:
%       - 0 : Both trees have the same structure.
%       - 1 : Different structures BUT obj structure is contained in 'otherTree' structure.
%       - 2 : Different structures BUT 'otherTree' structure is contained in obj structure.
%       - 3 : The trees have disjoint structures.
%
% - obj.genStructureTree():
%       Return a structure tree constructed from the tree obj.
%       A structure tree is defined so that every node's key is equal to
%       the number of leaves that descend from that node (e.g. the root
%       node's key is equal to the number of leaves of the tree, and each
%       leaf node has a key value of 1).
%
% - obj.structure2bitstream(bitstream):
%       Convert the tree structure to binary data and add it to the input
%       'bitstream' (BitStream object).
%       Only the tree structure is converted. The nodes' keys are ignored.s
%
% - bitstream2tree(bitstream):    (static function)
%       Returns a tree (with empty keys) constructed from a bitstream
%       containing tree structure (generated with the structure2bitstream
%       method). The input 'bitstream' is a BitStream object.
%

classdef BinNode < handle
    properties
        key
        left
        right
    end
    
    methods
        
        function obj = BinNode(key,left,right)
            obj.key = key;
            if(nargin>1)
                obj.left=left;
            end
            if(nargin>2)
                obj.right=right;
            end            
        end
        
        function b = isLeaf(obj)
            b = isempty(obj.left) && isempty(obj.right);
        end
        
        function n = numChildren(obj)
            n=0;
            if(~isempty(obj.left)), n=n+1;end
            if(~isempty(obj.right)), n=n+1;end
        end
        
        function d = depth(obj)
            d = obj.countLevels(0);
        end
        
        function n = numNodesAtLevel(obj,level)
            n = numNodesAtLevelRecur(obj,level,0,0);
        end
        
        function n = numLeavesUpToLevel(obj,level)
            n = numLeavesUpToLevelRecur(obj,level,0,0);
        end
        
        function n = numLeaves(obj)
            n = numLeavesRecur(obj,0);
        end
        
        function cList = listLeavesUpToLevel(obj,level)
            cList = listLeavesUpToLevelRecur(obj,level,0);
        end
        
        function tree = treeUpToLevel(obj,level)
            tree = treeUpToLevelRecur(obj,level,0);
        end
        
        function flag = structureDiff(obj,otherTree)
            flag = BinNode.structureDiffRecurr(obj,otherTree,0);
        end
        
        function structTree = genStructureTree(obj)
            if(~isempty(obj))
                if(obj.isLeaf())
                    structTree = BinNode(1);
                else
                    if(~isempty(obj.left))
                        NodeL = genStructureTree(obj.left);
                    else
                        NodeL = [];
                    end
                    if(~isempty(obj.right))
                        NodeR = genStructureTree(obj.right);
                    else
                        NodeR = [];
                    end
                    structTree = BinNode(NodeL.key + NodeR.key, NodeL, NodeR);
                end
            else
                structTree = [];
            end
        end
        
        function structure2bitstream(obj,bitstream)
            queue = obj;
            while ~isempty(queue)
                curNode = queue(1);
                queue(1)=[];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(curNode.isLeaf())
                    bitstream.add_bin(0);
                else
                    bitstream.add_bin(1);
                    if(curNode.numChildren==1), error('A tree node with single children was found. This is forbidden when converting tree structure to bitstream!');end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    queue(end+1) = curNode.left;
                    queue(end+1) = curNode.right;
                end
            end
        end
        
    end
    
    methods ( Static )
         function tree = bitstream2tree(bitstream)
            tree = BinNode([]);
            queue = tree;
            while ~bitstream.isempty() && ~isempty(queue)
                curNode = queue(1);
                queue(1)=[];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(bitstream.read_bin())
                    curNode.left = BinNode([]);
                    curNode.right = BinNode([]);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    queue(end+1) = curNode.left;
                    queue(end+1) = curNode.right;
                end
            end
            if(~isempty(queue))
                error('Bitstream incompatible with binary tree structure data (may be a truncated bitstream)');
            end
            
         end
    end
    
    methods (Access=private)
        
        function n = countLevels(obj,n)
            if(~isempty(obj))
                if(isempty(obj.left))
                    nLeft = n;
                else
                    nLeft = obj.left.countLevels(n+1);
                end
                if(isempty(obj.right))
                    nRight = n;
                else
                    nRight = obj.right.countLevels(n+1);
                end
                n = max( nLeft, nRight );
            end
        end
        
        function n = numNodesAtLevelRecur(obj,tgtLvl,curLvl,n)
            if(~isempty(obj))
                if(tgtLvl==curLvl)
                    n = n+1;
                else
                    if(~isempty(obj.left))
                        n = numNodesAtLevelRecur(obj.left,tgtLvl,curLvl+1,n);
                    end
                    if(~isempty(obj.right))
                        n = numNodesAtLevelRecur(obj.right,tgtLvl,curLvl+1,n);
                    end
                end
            end
        end
        
        function n = numLeavesUpToLevelRecur(obj,tgtLvl,curLvl,n)
            if(~isempty(obj))
                if(tgtLvl==curLvl || obj.isLeaf())
                    n = n+1;
                else
                    if(~isempty(obj.left))
                        n = numLeavesUpToLevelRecur(obj.left,tgtLvl,curLvl+1,n);
                    end
                    if(~isempty(obj.right))
                        n = numLeavesUpToLevelRecur(obj.right,tgtLvl,curLvl+1,n);
                    end
                end
            end
        end
        
        function n = numLeavesRecur(obj,n)
            if(~isempty(obj))
                if(obj.isLeaf)
                    n = n+1;
                else
                    if(~isempty(obj.left))
                        n = numLeavesRecur(obj.left,n);
                    end
                    if(~isempty(obj.right))
                        n = numLeavesRecur(obj.right,n);
                    end
                end
            end
        end
        
        function cList = listLeavesUpToLevelRecur(obj,tgtLvl,curLvl)
            if(~isempty(obj))
                if(tgtLvl==curLvl || obj.isLeaf())
                    cList = {obj.key};
                else
                    if(~isempty(obj.left))
                        cListL = listLeavesUpToLevelRecur(obj.left,tgtLvl,curLvl+1);
                    else
                        cListL = {};
                    end
                    if(~isempty(obj.right))
                        cListR = listLeavesUpToLevelRecur(obj.right,tgtLvl,curLvl+1);
                    else
                        cListR = {};
                    end
                    cList = [cListL, cListR];
                end
            else
                cList={};
            end
        end
        
        function tree = treeUpToLevelRecur(obj,tgtLvl,curLvl)
            if(~isempty(obj))
                if(tgtLvl==curLvl || obj.isLeaf())
                    tree = BinNode(obj.key);
                else
                    
                    if(~isempty(obj.left))
                        NodeL = treeUpToLevelRecur(obj.left,tgtLvl,curLvl+1);
                    else
                        NodeL = [];
                    end
                    if(~isempty(obj.right))
                        NodeR = treeUpToLevelRecur(obj.right,tgtLvl,curLvl+1);
                    else
                        NodeR = [];
                    end
                    tree = BinNode(obj.key, NodeL, NodeR);
                    tree.key = obj.key;
                end
            else
                tree = [];
            end
        end        
        
    end
    
    methods (Access=private, Static)
        function flag = structureDiffRecurr(tree1,tree2,flag)
            if(~isempty(tree1) && isempty(tree2))
                if(flag==0 || flag==1)
                    flag=flag+2;
                 %   previously: similar (flag=0) ==> now: tree2 included in tree1 (flag=2).
                 %OR previously: tree1 included in tree2 (flag=1) ==> now: disjoint (flag=3).
                end
                
            elseif(isempty(tree1) && ~isempty(tree2))
                if(flag==0 || flag==2)
                    flag=flag+1;
                 %   previously: similar (flag=0) ==> now: tree1 included in tree2 (flag=1).
                 %OR previously: tree2 included in tree1 (flag=2) ==> now: disjoint (flag=3).
                end
            elseif(~isempty(tree1) && ~isempty(tree2))
                flag = BinNode.structureDiffRecurr(tree1.left,tree2.left,flag);
                if(flag>=3),return;end %We already know the two trees are disjoint ==> no need to explore further.
                flag = BinNode.structureDiffRecurr(tree1.right,tree2.right,flag);
            end
        end
    end
    
end