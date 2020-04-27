#!/usr/bin/python3
import collections
import logging
from PyQt5.QtCore import QAbstractItemModel, QModelIndex, Qt


class ProjectItemModel(QAbstractItemModel):
    """
    Implementation of QAbstractItemModel. Had to implement in order to
    allow for drag and drop of my custom items (which have a second field of data
    stored within).
    Has two user roles: sequence role, and window role.
    -->Sequence role marks a dictionary of sequences that are stored within the node.\n
        For the sequence tree, that should only ever contain one sequence.
    -->Window role marks an alignment window that is built and saved. It is associated \n
        with a dictionary of all sequences found within that window (usually, it will be aligned).
    """


    def __init__(self, parent=None):
        super(ProjectItemModel, self).__init__()
        if not parent:
            self.rootItem = SeqNode("Default")
        else:
            self.rootItem = parent

    def appendRow(self, child):
        self.rootItem.addChild(child)

    def data(self, index, role=Qt.DisplayRole):
        try:
            item = index.internalPointer()
        except IndexError:
            return Qt.QVariant()
        if role == Qt.DisplayRole:
            return item
        if role == self.SequenceRole:
            return item.sequence()
        if role == self.WindowRole:
            return item.window()
        return Qt.QVariant()

    def index(self, row=None, column=None, parent=None):
        if not self.hasIndex(row, parent):
            return QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        childItem = parentItem.child(row)
        if childItem:
            return self.createIndex(row, childItem)
        else:
            return QModelIndex()

    def parent(self, index):
        if not index.isValid():
            return QModelIndex()

        childItem = index.internalPointer()
        parentItem = childItem.parent()

        if parentItem == self.rootItem:
            return QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)

    def rowCount(self, parent=None):
        if parent.column() > 0:
            return 0

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        return parentItem.getChildCount()

    def columnCount(self, parent=None):
        # necessary because I'm implementing an abstract class
        # but currently not used
        return 0


class TreeNode(object):
    def __init__(self, name, sequence=None, window=None, parent=None):
        self.nodeLogger = logging.getLogger("Node")
        self._name = name
        self._window = window
        self._parent = parent
        self._children = []
        self._seqs = {}
        if not sequence:
            self._seqs = {name, None}
            self.nodeLogger.debug("Empty node created")
        elif isinstance(sequence, collections.Mapping):
            self.nodeLogger.debug("Detected prebuilt sequence dictionary")
            self._seqs = sequence
        elif isinstance(sequence, str):
            self.nodeLogger.debug("Single sequence detected -- created single sequence node.")
            self._seqs = {name, sequence}
        else:
            self.nodeLogger.debug("Something went wrong with node creation!")

    # GET FUNCTIONS
    def name(self):
        return self._name

    def window(self):
        if self.window:
            return self.window

    def sequences(self):
        if len(list(self._seqs.keys())) > 1:
            return self._seqs
        else:
            return list(self._seqs.values())[0]

    def parent(self):
        return self._parent

    def children(self):
        return self._children

    def child(self, row):
        return self._children[row]

    def childCount(self):
        return len(self._children)

    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0

    # SET FUNCTIONS
    def setParent(self, parent):
        self._parent = parent

    def appendChild(self, node):
        self._children.append(node)


class TreeItemModel(QAbstractItemModel):
    SequenceRole = Qt.UserRole + 1
    WindowRole = Qt.UserRole + 2
    _roles = {SequenceRole: b"seqs", WindowRole: b"window"}

    def __init__(self, data, parent=None):
        super(TreeItemModel, self).__init__(parent)

        self.rootItem = TreeNode("Default")
        #self.setupModelData(data.split('\n'), self.rootItem)

    def columnCount(self, parent):
        if parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return self.rootItem.columnCount()

    def data(self, index, role):
        if not index.isValid():
            return None

        if role != Qt.DisplayRole:
            return None

        item = index.internalPointer()

        return item.data(index.column())

    def flags(self, index):
        if not index.isValid():
            return Qt.NoItemFlags

        return Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.rootItem.data(section)

        return None

    def index(self, row, column, parent):
        if not self.hasIndex(row, column, parent):
            return QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        childItem = parentItem.child(row)
        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QModelIndex()

    def parent(self, index):
        if not index.isValid():
            return QModelIndex()

        childItem = index.internalPointer()
        parentItem = childItem.parent()

        if parentItem == self.rootItem:
            return QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)

    def rowCount(self, parent):
        if parent.column() > 0:
            return 0

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        return parentItem.childCount()

    def setupModelData(self, lines, parent):
        parents = [parent]
        indentations = [0]

        number = 0

        while number < len(lines):
            position = 0
            while position < len(lines[number]):
                if lines[number][position] != ' ':
                    break
                position += 1

            lineData = lines[number][position:].trimmed()

            if lineData:
                # Read the column data from the rest of the line.
                columnData = [s for s in lineData.split('\t') if s]

                if position > indentations[-1]:
                    # The last child of the current parent is now the new
                    # parent unless the current parent has no children.

                    if parents[-1].childCount() > 0:
                        parents.append(parents[-1].child(parents[-1].childCount() - 1))
                        indentations.append(position)

                else:
                    while position < indentations[-1] and len(parents) > 0:
                        parents.pop()
                        indentations.pop()

                # Append a new item to the current parent's list of children.
                parents[-1].appendChild(TreeItem(columnData, parents[-1]))

            number += 1
