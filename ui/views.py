#!/usr/bin/python3

from PyQt5.QtWidgets import QWidget, QMdiSubWindow, QMdiArea, QTabBar, QTabWidget
from PyQt5.QtCore import Qt, pyqtSignal
from ui import alignment_ui
import textwrap as tw

class MDIArea(QMdiArea):
    def __init__(self):
        super(MDIArea, self).__init__()
        self.setViewMode(1)  # tabbed
        self.tabBar = self.findChild(QTabBar)
        self.setupTabBar()
        self.lastClosedTab = 0

    def resizeEvent(self, event):
        for sub in self.subWindowList():
            sub.resizeEvent(event)
        super(MDIArea, self).resizeEvent(event)

    def setupTabBar(self):
        #self.tabBar.setAutoHide(True)
        self.setTabsMovable(True)
        self.setTabsClosable(True)
        self.tabBar.tabCloseRequested.connect(self.closeTab)
        #self.tabBar.tabMoved.connect(self.moveTab)

    def closeTab(self, index):
        self.lastClosedTab = index
        window = self.childWindows[index]
        self.childWindows.remove(window)
        #self.removeSubWindow(window)
        #self.tabBar.removeTab(index)
        self.activeSubWindow().showMaximized()

    #def moveTab(self, orig, new):
    #    print("Tab moved")


    def addSubWindow(self, window, flags=Qt.WindowFlags()):
        #self.childWindows.append(window)
        super(MDIArea, self).addSubWindow(window, flags)
        for sub in self.subWindowList():
            sub.showMinimized()
        window.show()
        self.setActiveSubWindow(window)
        self.activeSubWindow().showMaximized()

    def setActiveSubWindow(self, window):
        titles = []
        for index in range(len(self.tabBar.children())):
            titles.append(self.tabBar.tabText(index))
        print(titles)
        print(self.subWindowList())
        #print(self.childWindows)
        if window.windowTitle() not in titles:
            print("adding tab")
            self.addSubWindow(window)
            self.tabBar.addTab(window.windowIcon(), window.windowTitle())
            self.activeSubWindow().showMaximized()
        else:
            super(MDIArea, self).setActiveSubWindow(window)
            self.activeSubWindow().showMaximized()


class MDISubWindow(QMdiSubWindow):
    """
    Have to subclass QMdiSubWindow because it doesn't automatically
    show the widget if I close the window, which is strange and annoying.
    """
    def __init__(self):
        super(MDISubWindow, self).__init__()
        self.setAttribute(Qt.WA_DeleteOnClose, False)

    def show(self):
        self.widget().show()
        super(MDISubWindow, self).show()

    def close(self):
        print("Closing subwindow")
        super(MDISubWindow, self).close()
        #super(MDISubWindow, self).hide()

    def closeEvent(self, event):
        print("Close event on subwindow!")
        self.mdiArea().removeSubWindow(self)
        self.close()



class AlignSubWindow(QWidget, alignment_ui.Ui_aliWindow):
    """
    Alignment subwindow UI. Takes in a dictionary of sequences that have been aligned and arranges them.
    """
    resized = pyqtSignal()

    def __init__(self, sequences):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self._seqs = sequences
        self.resized.connect(self.seqArrange)
        self.alignPane.verticalScrollBar().valueChanged.connect(self.namePane.verticalScrollBar().setValue)

    def resizeEvent(self, event):
        self.resized.emit()
        return super(AlignSubWindow, self).resizeEvent(event)

    def seqArrange(self):
        splitseqs = []
        prettynames = []
        prettyseqs = []
        maxname = 0
        nseqs = len(self._seqs.keys())
        width = (self.alignPane.size().width())
        charpx = self.alignPane.fontMetrics().averageCharWidth()
        nlines = 0
        for name, seq in self._seqs.items():
            maxchars = round(width / charpx) - 2 # for the scroll bar
            lines = tw.wrap(seq, maxchars)
            splitseqs.append([name, lines])
            if len(lines) > nlines:
                nlines = len(lines)
            if len(name) > maxname:
                maxname = len(name)
        # Set name window width to account for max name
        self.namePane.setMinimumWidth((maxname*charpx)+5)

        # Adjust the number of lines so that it accounts for the number of sequences
        # as well as a blank line (minus one for the last blank line)
        nlines = nlines*(nseqs+1)-1
        subline = 0
        seqid = 0
        # This loop creates a single alignment array by threading each of the individual lines
        # into each other.
        # TODO: Add TEXT FORMATTING to this section!! May need store as alternative text storage besides array.
        # TODO: For example, consider adding to the text window directly rather than as an array!
        for line in range(nlines):
            if seqid != nseqs:
                try:
                    prettynames.append(splitseqs[seqid][0])
                    prettyseqs.append(splitseqs[seqid][1][subline])
                    seqid += 1
                except IndexError:
                    prettynames.append(splitseqs[seqid][0])
                    prettyseqs.append("")
                    seqid += 1
            elif seqid == nseqs:
                prettyseqs.append("")
                prettynames.append("")
                seqid = 0
                subline += 1
        self.alignPane.setText("\n".join(prettyseqs))
        self.namePane.setAlignment(Qt.AlignRight)
        self.namePane.setText(prettynames[0])
        for line in prettynames[1:]:
            self.namePane.setAlignment(Qt.AlignRight)
            self.namePane.append(line)
        self.namePane.verticalScrollBar().setValue(self.alignPane.verticalScrollBar().value())

    def seqs(self):
        return self._seqs

