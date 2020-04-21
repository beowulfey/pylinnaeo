#!/usr/bin/python3

# This is the main class for a "workspace"
# Sets up an empty tree and populates it.

from PyQt5.QtGui import QStandardItem


# Each node of the tree comprises a sequence name and its raw sequence info.
class WorkNode(QStandardItem):
    def __init__(self, item=None, seq=None, parent=None):
        super(WorkNode, self).__init__(item)
        self.parentItem = parent
        self.itemSeq = seq
        self.childItems = []

    def appendChild(self, item):
        self.childItems.append(item)

    def getChild(self, row):
        return self.childItems[row]

    def getChildCount(self):
        return len(self.childItems)

    def getSeq(self):
        return self.itemSeq

    def getParent(self):
        return self.parentItem

    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0


#class WorkspaceModel(QStandardItemModel):
#    def __init__(self, data, parent=None):
#        super(WorkspaceModel, self).__init__(parent)
