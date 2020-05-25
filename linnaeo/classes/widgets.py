#!/usr/bin/python3
import logging
import sys

from PyQt5.QtCore import Qt, pyqtSignal, QRegularExpression, QSize, QPoint
from PyQt5.QtGui import QStandardItemModel, QFont, QFontDatabase, QColor, QSyntaxHighlighter, QTextCharFormat, \
    QTextCursor, QFontMetricsF, QTextDocument, QCursor
from PyQt5.QtWidgets import QWidget, QMdiSubWindow, QMdiArea, QTabBar, QTreeView, QSizePolicy, QAbstractItemView, \
    QDialog, QDialogButtonBox, QApplication, QTextEdit, QAbstractScrollArea, QToolTip
from PyQt5.uic.properties import QtCore, QtWidgets

from linnaeo.resources import linnaeo_rc
from linnaeo.classes import utilities
from linnaeo.ui import alignment_ui, quit_ui


class TreeView(QTreeView):
    generalClick = pyqtSignal()

    def __init__(self):
        super(TreeView, self).__init__()
        self.sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setSizePolicy(self.sizePolicy)
        self.setMinimumWidth(150)
        self.setEditTriggers(QAbstractItemView.SelectedClicked)
        self.setDragEnabled(True)
        self.setDragDropMode(QAbstractItemView.InternalMove)

    def mousePressEvent(self, event):
        self.generalClick.emit()
        return super(TreeView, self).mousePressEvent(event)


class MDIArea(QMdiArea):
    """
    Custom QMdiArea class -- the native class had some strange bugs with closing tabs, where it would keep the tab
    due to how the closing was happening. I've subclassed it to make things easy in finding bugs.
    This works well and should be more customizable if needed.
    """
    def __init__(self):
        super(MDIArea, self).__init__()
        self.tabbed = False  # tabbed by default
        self.tabBar = None
        self.setTabs(True) if self.tabbed else self.setTabs(False)
        self.WindowOrder(2)

    def toggleTabs(self):
        self.tabbed = not self.tabbed
        if self.tabbed:
            self.setTabs(True)
        else:
            self.setTabs(False)
            for sub in self.subWindowList():
                sub.showNormal()
            self.tileSubWindows()

    def setTabs(self, on):
        if on:
            self.setViewMode(1)
            self.tabbed = True
            self.tabBar = self.findChild(QTabBar)
            self.setupTabBar()
        else:
            self.tabbed = False
            self.setViewMode(0)

    def setupTabBar(self):
        # self.tabBar.setAutoHide(True)
        self.setTabsMovable(True)
        self.setTabsClosable(True)
        self.tabBar.tabCloseRequested.connect(self.closeTab)

    def closeTab(self):
        try:
            self.activeSubWindow().showMaximized()
        except:
            # So I can close all tabs
            pass

    def resizeEvent(self, event):
        # TODO: CAN I DELETE THIS??
        # passes a resize event to all subwindows to make sure the sequence is updated
        for sub in self.subWindowList():
            sub.resizeEvent(event)
        return super(MDIArea, self).resizeEvent(event)

    def addSubWindow(self, window, flags=Qt.WindowFlags()):
        super(MDIArea, self).addSubWindow(window, flags)
        if self.tabbed:
            print(self.tabbed)
            for sub in self.subWindowList():
                sub.showMinimized()
            self.activeSubWindow().showMaximized()
        window.show()
        self.setActiveSubWindow(window)

    def setActiveSubWindow(self, window):
        if self.tabbed:
            titles = []
            for index in range(len(self.tabBar.children())):
                titles.append(self.tabBar.tabText(index))
            if window.windowTitle() not in titles:
                self.addSubWindow(window)
                self.tabBar.addTab(window.windowIcon(), window.windowTitle())
                self.activeSubWindow().showMaximized()
            else:
                super(MDIArea, self).setActiveSubWindow(window)
                self.activeSubWindow().showMaximized()
        else:
            super(MDIArea, self).setActiveSubWindow(window)
            window.widget().resized.emit()
            #print("Unable to resize")
            #self.activeSubWindow().showMaximized()


class MDISubWindow(QMdiSubWindow):
    """
    Have to subclass QMdiSubWindow because it doesn't automatically
    show the widget if I close the window, which is strange and annoying.
    """
    def __init__(self):
        super(MDISubWindow, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose, False)
        self._widget = None
        self.setMouseTracking(True)
        #self.installEventFilter(self.widget().alignPane)

        # TODO: TAKE OUT EXTRA CLOSE COMMAND IN MDISUBWINDOW
        # remove extra close command
        # menu = self.systemMenu()
        # for action in menu.actions():
        #    if action.text()=="&Close":
        #        print("Found")
        #        menu.actions().remove(action)
        # self.setSystemMenu(menu)

    #def eventFilter(self, obj, event):
    #    if event.type() == 129:
    #        print("HOVER")
    #    if event.type() == 5:
    #        print("MOUSE MOVE")
    #    return super().eventFilter(obj, event)

    def mouseMoveEvent(self, event):
        #print(QCursor.pos())
        return super().mouseMoveEvent(event)

    def event(self, event):
        # EventFilter doesn't capture type 2 events on title bar of subwindow for some reason
        # Is there a better way to do this???
        #print(event, event.type())
        if event.type() == 2:
            #linnaeo = self.parentWidget().parentWidget().parentWidget().parentWidget().parentWidget().parentWidget()
            #print("REDRAWING FROM MDI")
            self._widget.userIsResizing = True
            #self._widget.seqArrangeNoColor()
        elif event.type() == 3:
            #print("DONE REDRAWING FROM MDI")
            self._widget.userIsResizing = False
            self._widget.resized.emit()
            #self._widget.seqArrangeColor()
        return super().event(event)

    def setWidget(self, widget):
        self._widget = widget
        super(MDISubWindow, self).setWidget(widget)

    def widget(self):
        return self._widget

    def show(self):
        self._widget.show()
        super(MDISubWindow, self).show()

    def closeEvent(self, event):
        if self.mdiArea():
            self.mdiArea().removeSubWindow(self)
        self.close()
        return super(MDISubWindow, self).closeEvent(event)

    def close(self):
        super(MDISubWindow, self).close()


class ItemModel(QStandardItemModel):
    nameChanging = pyqtSignal()
    dupeName = pyqtSignal()
    nameChanged = pyqtSignal()

    def __init__(self, windows, seqTree=False):
        super(ItemModel, self).__init__()
        self.modelLogger = logging.getLogger("ItemModel")
        self.lastClickedNode = None
        #self._root = root
        self._windows = windows
        self._titles = []
        self.isSeqs = False
        if seqTree:
            self.isSeqs = True

    def updateLastClicked(self, node):
        self.lastClickedNode = node

    def updateNames(self, titles):
        self._titles = titles

    def updateWindows(self, windows):
        self._windows = windows

    def setData(self, index, value, role=Qt.UserRole+1):
        oldvalue = None
        self.modelLogger.debug("Updating data for node")
        if self.isSeqs:
            self.modelLogger.debug("Sequence node; checking name!")
            #self.nameChanging.emit()
            # Only do this check if this is coming from the top Tree and is not a folder
            if value != self.lastClickedNode.text():
                oldvalue = self.lastClickedNode.text()
                newvalue, self._titles = utilities.checkName(value, self._titles)
                if newvalue != value:
                    self.modelLogger.debug("Item duplicates a different node! "+str(value)+" in "+str(self._titles))
                    value = newvalue
                    self.dupeName.emit()
                    self.modelLogger.debug("Name changed to "+str(value))
                try:
                    print("setDataTryLoop")
                    sub = self._windows[self.itemFromIndex(index).data(role=Qt.UserRole+3)]
                    sub.widget().nameChange.emit(oldvalue, newvalue)
                except KeyError:
                    print("KEY ERROR")
                    pass
                if self.itemFromIndex(index).data(role=Qt.UserRole+2):
                    seqr = self.itemFromIndex(index).data(role=Qt.UserRole+2)[0]
                    seqr.name = newvalue
        self.modelLogger.debug("Setting node text to "+str(value))
        self.itemFromIndex(index).setData(value)
        self.itemFromIndex(index).setText(value)
        self.nameChanged.emit()
        try:
            sub = self._windows[self.itemFromIndex(index).data(role=Qt.UserRole+3)]
            sub.setWindowTitle(value)

            sub.mdiArea().setActiveSubWindow(sub)
        except:
            exctype, val = sys.exc_info()[:2]
            self.modelLogger.debug("Detected exception, but probably fine: "+str(exctype))
        finally:
            return super(ItemModel, self).setData(index, value, role)


class AlignDoc(QTextDocument):
    def __init__(self):
        super().__init__()
        self.setMouseTracking(True)

    def event(self, event):
        print(event, event.type())
        return super().event(event)

    def mouseMoveEvent(self, event):
        print(QCursor.pos())
        return super().mouseMoveEvent(event)


class AlignPane(QTextEdit):

    toolTipReq = pyqtSignal(QPoint, str)

    def __init__(self, parent):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.tracking = False
        self.offset = (334,208)

        sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.toolTipReq.connect(self.getSeqTT)
        self.cursorPositionChanged.connect(self.getPos)

        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(200, 100))
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)
        self.setLineWrapMode(QTextEdit.NoWrap)
        self.setReadOnly(True)
        self.setTextInteractionFlags(Qt.TextSelectableByKeyboard | Qt.TextSelectableByMouse)
        self.setObjectName("alignPane")
        self.setCursorWidth(1)
        self.setToolTipDuration(-1)

    #def event(self, event):
        #print(event, event.type())
    #    if event.type() == 110:
    #        self.toolTipReq.emit()
    #    return super().event(event)

    def getPos(self):
        pass
        #print("POS CHANGED")
        #print(self.textCursor().position())
        #print(self.textCursor().positionInBlock())

    def getSeqTT(self, pos, selected):
        print("TOOLTIP REQ")
        #self.setToolTip()
        QToolTip.showText(pos, selected)
        #self.unsetCursor()
        self.textCursor().clearSelection()

    #def toolTipEvent(self, event):
    #    print("TOOLTIP")
    #    self.setToolTip(self.textCursor().selectedText())
    #    self.textCursor().clearSelection()
    #    return super().toolTipEvent(event)

    """def mouseMoveEvent(self, event):
        if self.tracking:
            self.unsetCursor()
            self.cursorForPosition(event.globalPos())
            self.mousePressEvent(event)
        else:
            super().mouseMoveEvent(event)


    def mousePressEvent(self, event):
        self.unsetCursor()
        super().mousePressEvent(event)
        print("CLICK!", event.globalPos())
        self.setTextCursor(self.cursorForPosition(self.mapFromGlobal(event.globalPos())))
        self.moveCursor(QTextCursor.Left, mode=QTextCursor.MoveAnchor)
        self.moveCursor(QTextCursor.Left, mode=QTextCursor.MoveAnchor)
        self.moveCursor(QTextCursor.Left, mode=QTextCursor.KeepAnchor)
        text = self.textCursor().selectedText()
        self.unsetCursor()

        self.toolTipReq.emit(event.globalPos(), text)
        self.tracking = True

    def mouseReleaseEvent(self, event):
        print("UNCLICK")
        super().mouseReleaseEvent(event)
        self.tracking = False





    #def mouseMoveEvent(self, event):
        # TODO THIS DOES NOT MATCH TRUE CURSOR ... QCURSOR RETURNS ABSOLUTE VALUE!
        #cursor = self.mapFromGlobal(QPoint(QCursor.pos()))
        #tcursor = self.cursorForPosition(cursor)
        #self.setTextCursor(tcursor)
        #print(QCursor.pos())
        #return super().mouseMoveEvent(event)"""



