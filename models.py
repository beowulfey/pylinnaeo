#!/usr/bin/python3

# This is the main class for a "workspace"
# Sets up an empty tree and populates it.

from PyQt5.QtGui import QStandardItem


# Each node of the tree comprises a sequence name and its raw sequence info.
class SeqNode(QStandardItem):
    def __init__(self, item=None, seq=None, parent=None):
        super(SeqNode, self).__init__(item)
        self.parentItem = parent
        self.itemSeq = seq
        self.childItems = []

    def getSeq(self):
        return self.itemSeq

    """def appendChild(self, item):
        self.childItems.append(item)

    def getChild(self, row):
        return self.childItems[row]

    def getChildCount(self):
        return len(self.childItems)

    def getParent(self):
        return self.parentItem

    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0"""

# Not sure if it's necessary to have a separate class for the workspace.
# They functionally work the same; this one just stores a window instead of a sequence.
# May be good for detecting different types later?...

class WorkspaceNode(SeqNode):
    def __init__(self, item=None, window=None, parent=None):
        super(WorkspaceNode, self).__init__(item)
        self.parentItem = parent
        self.window = window
        self.childItems = []

    def getWindow(self):
        return self.window
