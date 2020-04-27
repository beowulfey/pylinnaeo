#!/usr/bin/python3

from PyQt5.QtWidgets import QWidget, QMdiSubWindow
from ui import alignment_ui


class AlignSubWindow(QWidget, alignment_ui.Ui_Form):
    def __init__(self, sequences):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self._seqs = sequences
        for seq in self._seqs:
            self.textEdit.setText(str(self._seqs[seq]))
            self.textEdit_2.setText(seq)

    def getSeqs(self):
        return self._seqs


class MDISubWindow(QMdiSubWindow):
    def __init__(self):
        super(MDISubWindow, self).__init__()

    def show(self):
        self.widget().show()
        super(MDISubWindow, self).show()
