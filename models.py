#!/usr/bin/python3

# This is the main class for a "workspace"
# Sets up an empty tree and populates it.


from PyQt5.QtGui import QStandardItem
from PyQt5.QtCore import QAbstractItemModel, QModelIndex, Qt


# Each node of the tree comprises a sequence name and its raw sequence info.
class SeqNode(QStandardItem):
    """
    Extension of QStandardItem that contains additional data
    Functionally identical to WorkspaceNode below; may merge.
    """

    def __init__(self, name=None, sequence=None, window=None, parent=None):
        super(SeqNode, self).__init__(name)
        self.FolderRole = Qt.UserRole + 1
        self.SequenceRole = Qt.UserRole + 2
        self.WindowRole = Qt.UserRole + 3
        # Set Default Flags
        self.setFlags(Qt.ItemIsEnabled | Qt.ItemIsEditable | Qt.ItemIsSelectable |
                      Qt.ItemIsUserCheckable | Qt.ItemIsDragEnabled)
        #self._name = name
        #self._parent = parent
        #self._childItems = []
        if sequence:
            super(SeqNode, self).setData(sequence, self.SequenceRole)
        if window:
            #self.setFlags(self.flags() or self.setFlags(Qt.ItemIsDropEnabled))
            super(SeqNode, self).setData(window, self.WindowRole)

    def clone(self):
        return super(SeqNode, self).clone()

    def sequence(self):
        return super(SeqNode, self).data(self.SequenceRole)

    def window(self):
        return super(SeqNode, self).data(self.WindowRole)

    def child(self, row, column=None):
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
        return 0

# Not sure if it's necessary to have a separate class for the workspace.
# They functionally work the same; this one just stores a window instead of a sequence.
# May be good for detecting different types later?...
