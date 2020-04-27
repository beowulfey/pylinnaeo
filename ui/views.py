#!/usr/bin/python3

from PyQt5.QtWidgets import QWidget, QMdiSubWindow
from PyQt5.QtCore import Qt, pyqtSignal
from ui import alignment_ui


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


class AlignSubWindow(QWidget, alignment_ui.Ui_Form):
    """
    Alignment subwindow UI. Takes in a dictionary of sequences that have been aligned and arranges them.
    """
    resized = pyqtSignal()

    def __init__(self, sequences):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self._seqs = sequences
        for seq in self._seqs:
            self.alignPane.setText(str(self._seqs[seq]))
            self.namePane.setText(seq)
        self.namePane.setAlignment(Qt.AlignRight)
        self.resized.connect(self.seqArrange)
        #self.alignWidget.resizeEvent.connect(self.seqArrange)

    def resizeEvent(self, event):
        self.resized.emit()
        return super(AlignSubWindow, self).resizeEvent(event)

    def seqArrange(self):
        print("RESIZE DETECTED!")
        print("Current width: "+ str(self.alignWidget.width()))

    def seqs(self):
        return self._seqs
