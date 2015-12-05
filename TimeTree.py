#!/usr/bin/env python
"""
    TimeTree
    --------

    Basic rooted phylogenetic time tree manipulation module which includes Tree
    and Node classes.  The Tree class can be initialized from a Newick string
    and provides methods for producing visualizations using either ASCII art or
    Matplotlib. (The latter requires matplotlib and scipy.)

    :copyright: 2015 by Tim Vaughan
    :license: GPL version 3.0, see http://www.gnu.org/licenses/gpl.html for more details
"""

import re

class Node:
    """A node in a phylogenetic tree."""

    def __init__(self):
        self.height = None
        self.time = None
        self.branchLength = None
        self.parent = None
        self.children = []
        self.annotations = {}
        self.label = None

    def isRoot(self):
        """Returns true if this node is the root of a tree."""

        return self.parent == None

    def isLeaf(self):
        """Returns true if this node is a leaf of a tree."""

        return len(self.children)==0

    def addChild(self, newChild):
        self.children.append(newChild)
        newChild.parent = self

    def getAllChildren(self):
        """Retrieve a list including this node and all of its decendents."""

        childList = [self]
        for child in self.children:
            childList.extend(child.getAllChildren())

        return childList

    def getLeaves(self):
        """Retrieve the list of leaves descending from this node."""

        if self.isLeaf():
            return [self]
        else:
            leaves = []
            for child in self.children:
                leaves.extend(child.getLeaves())
            return leaves

    def getNewick(self):
        """Retrieve a Newick representation of the tree below this node."""

        newick = ""
        if len(self.children)>0:
            newick += "("
            for i,v in enumerate(self.children):
                if i>0:
                    newick += ","
                newick += v.getNewick()
            newick += ")"

        if self.label != None:
            newick += self.label

        if len(self.annotations)>0:
            newick += '[&'
            isFirst = True
            for k,v in self.annotations.items():
                if isFirst:
                    isFirst = False
                else:
                    newick += ","
                newick += '{}="{}"'.format(k,v)
            newick += ']'

        if self.parent == None:
            if self.origin != None:
                newick += ":{}".format(self.origin - self.height)
            else:
                newick += ":0.0"
        else:
            newick += ":{}".format(self.parent.height - self.height)

        return newick

    def computeTimes(self, offset):
        """Set up time attribute of this node and each node below it.  The offset
        parameter indicates the time of this node's parent."""

        self.time = offset + self.branchLength
        for child in self.children:
            child.computeTimes(self.time)


class Tree:
    """A phylogenetic tree."""

    def __init__(self, arg):
        """Create a new phylogenetic tree.  The argument can either be a root
        node, a newick-formatted string, or a file object corresponding to an
        open newick or nexus-formatted file."""

        if type(arg) is Node:
            self.root = arg
            return

        # Need to parse string or string from file
        newickString = ""
        if type(arg) is str:
            newickString = arg.strip()

        if type(arg) is file:
            firstLine = arg.readline()

            newickString = ""
            if firstLine.lower().startswith('#nexus'):
                for line in arg.readlines():
                    if line.strip().lower().startswith("tree "):
                        newickString = line[(line.find("=")+1):].strip()
                        break
            else:
                newickString = firstLine.strip()

        self.root = self.loadFromString(newickString)


    def __repr__(self):
        return "Phylogenetic tree with {} nodes (including {} leaves).".format(len(self.getNodes()), len(self.getLeaves()))


    # Basic tree queries

    def getNodes(self):
        """Constructs and returns a list of all nodes in the tree."""

        return self.root.getAllChildren()

    def getLeaves(self):
        """Constructs and returns a list of all children in the tree."""

        return self.root.getLeaves()

    # Tree manipulation

    def sort(self, increasing=True):
        """Sort tree (in place) by ordering children according to the size of their subtrees."""

        def sortSubTree(node):
            thisCladeSize = 1
            for child in node.children:
                child._cladeSize = sortSubTree(child)
                thisCladeSize += child._cladeSize

            node.children.sort(key=lambda child: child._cladeSize, reverse=not increasing)

            return thisCladeSize

        sortSubTree(self.root)

        return self

    # Tree parsing code

    class ParseError(Exception):
        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    class ParseContext:
        def __init__(self, tokenList, valueList):
            self.tokenList = tokenList
            self.valueList = valueList
            self.idx = 0

        def acceptToken(self, token, manditory=False):
            if self.tokenList[self.idx] == token:
                self.idx += 1
                return True
            else:
                if not manditory:
                    return False
                else:
                    raise ParseError("Error parsing token {} ({})".format(self.tokenList[self.idx], self.valueList[self.idx]))

        def getLastValue(self):
            return self.valueList[self.idx-1]


    def loadFromString(self, string):
        """Parse string representation of tree and store in this object.

        Note: This is a utility method called from __init__()."""

        # Lexical analysis

        tokens = [
                ('LPAREN',  '\('),
                ('RPAREN',  '\)'),
                ('COLON',   ':'),
                ('STRING', '"[^"]*"'),
                ('STRING', '\'[^\']*\''),
                ('STRING',   '[a-zA-Z0-9_.-]+'),
                ('OPENA', '\[&'),
                ('EQUALS', '='),
                ('CLOSEA', '\]'),
                ('COMMA',   ','),
                ('SEMI',    ';')
                ]

        idx = 0
        tokenList = []
        valueList = []

        while idx < len(string):

            noMatch = True

            for k in range(len(tokens)):
                match = re.match(tokens[k][1], string[idx:])

                if match != None:
                    tokenList.append(tokens[k][0])
                    idx += len(match.group(0))

                    if tokens[k][0] == 'STRING':
                        value = match.group(0)
                        if len(value)>2:
                            if value[0]=='"' and value[len(value)-1]=='"':
                                value = value[1:(len(value)-1)]
                            else:
                                if value[0]=="'" and value[len(value)-1]=="'":
                                    value = value[1:(len(value)-1)]
                        valueList.append(value)
                    else:
                        valueList.append(None)

                    noMatch = False
                    break

            if noMatch:
                raise Tree.ParseError('Unrecognized character at position ' + str(idx) + ': \'' + string[idx] + '\'')

        # Recursive decent parser

        def ruleN(parent, ctx, depth):
            node = Node()
            if parent != None:
                parent.addChild(node)

            ruleS(node, ctx, depth)
            ruleL(node, ctx, depth)
            ruleA(node, ctx, depth)
            ruleB(node, ctx, depth)

            #print " "*depth + str(node.label) + ":" + str(node.branchLength)
            return node

        def ruleS(node, ctx, depth):
            if ctx.acceptToken('LPAREN'):
                #print " "*depth + "("
                ruleN(node, ctx, depth+1)
                ruleQ(node, ctx, depth)
                ctx.acceptToken('RPAREN', manditory=True)
                #print " "*depth + ")"

        def ruleQ(node, ctx, depth):
            if ctx.acceptToken('COMMA'):
                ruleN(node, ctx, depth+1)
                ruleQ(node, ctx, depth)

        def ruleL(node, ctx, depth):
            if ctx.acceptToken('STRING'):
                node.label = ctx.getLastValue()

        def ruleA(node, ctx, depth):
            if ctx.acceptToken('OPENA'):
                ruleC(node, ctx, depth)
                ruleD(node, ctx, depth)
                ctx.acceptToken('CLOSEA', manditory=True)

        def ruleC(node, ctx, depth):
            ctx.acceptToken('STRING', manditory=True)
            key = ctx.getLastValue()
            ctx.acceptToken('EQUALS', manditory=True)
            ctx.acceptToken('STRING', manditory=True)
            value = ctx.getLastValue()

            node.annotations[key] = value

        def ruleD(node, ctx, depth):
            if ctx.acceptToken('COMMA'):
                ruleC(node, ctx, depth)
                ruleD(node, ctx, depth)

        def ruleB(node, ctx, depth):
            if ctx.acceptToken('COLON'):
                ctx.acceptToken('STRING', manditory=True)
                node.branchLength = float(ctx.getLastValue())
            else:
                node.branchLength = 1.0

        ctx = Tree.ParseContext(tokenList, valueList)
        root = ruleN(None, ctx, 0)
        ctx.acceptToken("SEMI", manditory=True)

        # Compute node heights and times:

        root.computeTimes(0.0)

        maxTime = 0.0
        for node in root.getAllChildren():
            maxTime = max(maxTime, node.time)

        for node in root.getAllChildren():
            node.height = maxTime - node.time

        root.origin = root.height + root.branchLength

        return root


    # Visualization code

    def plot(self, **kwargs):
        """Simple tree visualization.  (Requires matplotlib.)"""

        from matplotlib.pylab import plot, text, gca, xlim, ylim, show

        if "color" not in kwargs.keys():
            kwargs['color'] = 'black'

        leaves = self.getLeaves()
        nodes = self.getNodes()
        pos = [0.]*len(nodes)

        def computePos(node):
            idx = nodes.index(node)
            if node.isLeaf():
                pos[idx] = leaves.index(node)
            else:
                for child in node.children:
                    pos[idx] += computePos(child)
                pos[idx] /= len(node.children)

            return pos[idx]
        computePos(self.root)

        for i in range(len(nodes)):

            if nodes[i].parent == None:
                plot([nodes[i].height, nodes[i].height+nodes[i].branchLength],
                        [pos[i], pos[i]], **kwargs)
            else:
                pi = nodes.index(nodes[i].parent)
                plot([nodes[i].height, nodes[i].parent.height, nodes[i].parent.height],
                        [pos[i], pos[i], pos[pi]], **kwargs)

            if nodes[i].label != None:
                text(nodes[i].height, pos[i], nodes[i].label, {'ha':'left', 'va':'bottom'}, rotation=45)

        for i in range(len(nodes)):
            if len(nodes[i].children)==1:
                plot([nodes[i].height], [pos[i]], 'ro')

        ylim(-0.05*len(leaves), 1.05*len(leaves)-1)

        ax = gca()
        ax.invert_xaxis()
        ax.yaxis.set_visible(False)
        ax.set_frame_on(False)
        ax.grid()

    def plot_ascii(self, width=70, labelLeaves=True):
        """Display crude ASCII representation of tree."""

        leaves = self.getLeaves()
        nodes = self.getNodes()
        pos = [0.]*len(nodes)

        grid = [[" " for i in range(width)] for leaf in leaves]

        def computePos(node):
            idx = nodes.index(node)
            if node.isLeaf():
                pos[idx] = leaves.index(node)
            else:
                for child in node.children:
                    pos[idx] += computePos(child)
                pos[idx] /= len(node.children)

            return pos[idx]
        computePos(self.root)

        # Edges
        for i in range(len(nodes)):

            if nodes[i].parent == None:
                x1 = int(nodes[i].height/self.root.origin*(width-1))
                x2 = width-1
                y = int(pos[i])
                grid[y][x1:x2] = ['-']*(x2-x1)
            else:
                pi = nodes.index(nodes[i].parent)
                x1 = int(nodes[i].height/self.root.origin*width)
                x2 = int(nodes[i].parent.height/self.root.origin*(width-1))
                y1 = int(pos[i])
                y2 = int(pos[pi])
                ymin = min(y1, y2)
                ymax = max(y1, y2)
                grid[y1][x1:x2] = ['-']*(x2-x1)
                for y in range(ymin+1, ymax):
                    grid[y][x2] = '|'

        # Corners
        for i in range(len(nodes)):

            if nodes[i].parent != None:
                pi = nodes.index(nodes[i].parent)
                x1 = int(nodes[i].height/self.root.origin*(width-1))
                x2 = int(nodes[i].parent.height/self.root.origin*(width-1))
                y1 = int(pos[i])
                y2 = int(pos[pi])
                ymin = min(y1, y2)
                ymax = max(y1, y2)
                grid[y1][x1:x2] = ['-']*(x2-x1)
                grid[ymin][x2] = "\\"
                grid[ymax][x2] = "/"

        # Nodes
        for i in range(len(nodes)):
            x = int(nodes[i].height/self.root.origin*(width-1))
            y = int(pos[i])
            if nodes[i].isLeaf():
                grid[y][x] = '*'
            else:
                grid[y][x] = '+'

        grid.reverse()
        leaves.reverse()
        for i in range(len(leaves)):
            row = grid[i]
            row.reverse()
            rowStr = "".join(row).rstrip()

            if labelLeaves and leaves[i].label != None:
                rowStr += " " + leaves[i].label

            print rowStr
