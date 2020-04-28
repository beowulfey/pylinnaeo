#!/usr/bin/python3

from PyQt5.QtWidgets import QWidget, QMdiSubWindow
from PyQt5.QtCore import Qt, pyqtSignal
from ui import alignment_ui
from PIL import ImageFont
import textwrap as tw


class MDISubWindow(QMdiSubWindow):
    """
    Have to subclass QMdiSubWindow because it doesn't automatically
    show the widget if I close the window, which is strange and annoying.
    """
    def __init__(self):
        super(MDISubWindow, self).__init__()

    def show(self):
        self.widget().show()
        super(MDISubWindow, self).show()


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
        #self.alignWidget.resizeEvent.connect(self.seqArrange)

    def resizeEvent(self, event):
        self.resized.emit()
        return super(AlignSubWindow, self).resizeEvent(event)

    def seqArrange(self):
        splitseqs = []
        prettynames = []
        prettyseqs = []
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

        # Adjust the number of lines so that it accounts for the number of sequences
        # as well as a blank line (minus one for the last blank line)
        nlines = nlines*(nseqs+1)-1
        subline = 0
        seqid = 0
        for line in range(nlines):
            #print(line, subline, seqid)
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
        self.namePane.setText("\n".join(prettynames))










        #self.namePane.setText(name)
        #self.namePane.setAlignment(Qt.AlignRight)
        #self.alignPane.setText(firstline)







    def seqs(self):
        return self._seqs
