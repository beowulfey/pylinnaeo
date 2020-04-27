#!/usr/bin/python3

# This is the main class for a "workspace"
# Sets up an empty tree and populates it.


from PyQt5.QtGui import QStandardItem
from PyQt5.QtCore import QAbstractItemModel, QModelIndex, Qt


# Each node of the tree comprises a sequence name and its raw sequence info.


class SeqNode(QStandardItem):
    """
    Extension of QStandardItem that contains additional data
    """
    def __init__(self, name=None, sequence=None, window=None, parent=None):
        super(SeqNode, self).__init__(name)
        self.SequenceRole = Qt.UserRole + 1
        self.WindowRole = Qt.UserRole + 2
        if sequence:
            print(super(SeqNode, self).flags())
            self.setFlags(self.flags() ^ Qt.ItemIsDropEnabled)
            super(SeqNode, self).setData(sequence, self.SequenceRole)
        if window:
            self.setFlags(self.flags() ^ Qt.ItemIsDropEnabled)
            super(SeqNode, self).setData(window, self.WindowRole)

    def sequence(self):
        return super(SeqNode, self).data(self.SequenceRole)

    def window(self):
        return super(SeqNode, self).data(self.WindowRole)

    def clone(self):
        return super(SeqNode, self).clone()

    """def child(self, row, column=None):
        return self._childItems[row]

    def addChild(self, item):
        self._childItems.append(item)

    def childCount(self):
        return len(self._childItems)

    def parent(self):
        return self._parent

    def row(self):
        if self._parent:
            return self._parent.childItems.index(self)
        return 0"""

# Not sure if it's necessary to have a separate class for the workspace.
# They functionally work the same; this one just stores a window instead of a sequence.
# May be good for detecting different types later?...
