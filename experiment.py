from collections import deque

class TreeNode:
    def __init__(self, value=0, left=None, right=None, parent=None):
        self.value = value
        self.left = left
        self.right = right
        self.parent = parent

def build_tree_from_preorder(preorder):
    root = TreeNode()
    prevbto = False
    prevbtc = False
    comma = False
    for x in preorder:
        if (x.isnumeric()):
            root.value = x
            if (comma):
                root = root.parent
                comma = False
        elif (x == '('):
            temp = TreeNode()
            temp.parent = root
            root.left = temp
            root.left = root
            prevbto = True
        elif (x == ','):
            temp = TreeNode()
            root = root.parent
            temp.parent = root
            root.right = temp
            root = root.right
            comma = True
        elif (x == ')'):
            root = root.parent
            prevbtc = True



# Example usage:
if __name__ == "__main__":
    # Pre-order traversal list including 'None' to mark the absence of a node
    preorder = '3(2(6,8),3(1,4))'

    ret = []
    for x in preorder:
        ret.append(x)
    print(ret)
    
    # Build the tree from pre-order traversal
    root = build_tree_from_preorder(ret)
    print(root)
    print("Level Order Traversal:")
    level_order_traversal(root)
