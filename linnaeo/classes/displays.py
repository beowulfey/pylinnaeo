from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor, QFontDatabase, QFont, QFontMetricsF, QTextCursor
from PyQt5.QtWidgets import QWidget, QApplication, QDialog, QDialogButtonBox

from linnaeo.resources import linnaeo_rc
from linnaeo.ui import alignment_ui, quit_ui


class AlignSubWindow(QWidget, alignment_ui.Ui_aliWindow):
    """
    Alignment SubWindow UI. Takes in a dictionary of sequences that have been aligned and arranges them.
    # TODO: This does not maintain the alignment order. Shant be helped?...
    """
    resized = pyqtSignal()


    # THEMES #

    pos = QColor(100, 140, 255)
    neg = QColor(255, 70, 90)
    cys = QColor(255, 255, 85)
    aro = QColor(145, 255, 168)

    defTheme = {
        # Charged; positive
        "R": pos, "H": aro, "K": pos,
        # Charged, negative
        "D": neg, "E": neg,
        # Misc
        "C": cys, #"G": gly, "A": ala,
        # Aromatic
        "W": aro, "F": aro, "Y": aro
    }

    def __init__(self, seqs):
        super(self.__class__, self).__init__()
        self.userIsResizing = False
        self.setupUi(self)
        self._seqs = seqs
        self.resized.connect(self.externalResizeDone)
        self.alignPane.verticalScrollBar().valueChanged.connect(self.namePane.verticalScrollBar().setValue)
        self.theme = None
        self.splitNames = []
        self.splitSeqs = []

        # FANCY FONTWORK
        # This maintains the font within the application.
        fid = QFontDatabase.addApplicationFont(":/fonts/LiberationMono.ttf")
        family = QFontDatabase.applicationFontFamilies(fid)[0]
        font = QFont(family, 10)
        self.fmF = QFontMetricsF(font)  # FontMetrics Float... because default FontMetrics gives Int
        self.alignPane.setFont(font)
        self.alignPane.setStyleSheet("QTextEdit {padding-left:20px; padding-right:0px; padding-top:0px; background-color: \
                                     rgb(255,255,255)}")
        self.namePane.setStyleSheet("QTextEdit {padding-top:0px;}")
        self.alignPane.setCursorWidth(0)

        self.refseq = None
        self.maxlen = 0


        # options to do
        # TODO: Implement these
        self.showRuler = False
        self.showColors = True
        self.relColors = False

        if self.showColors:
            self.theme = self.defTheme

        self.seqInit()

    def setTheme(self, theme):
        self.theme = theme

    def toggleRulers(self):
        self.showRuler = not self.showRuler

    def resizeEvent(self, event):
        """
        This gets called anytime the window is in the process of being redrawn. If the MDI Subwindow is maximized,
        it calls a resizeEvent upon release too.
        """
        if self.userIsResizing:
            #print("REDRAW FROM INSIDE")
            self.seqArrangeNoColor()
        elif not self.userIsResizing:
            self.resized.emit()
            #print("DONE FROM INSIDE")
        super(AlignSubWindow, self).resizeEvent(event)
        # self.oldwidth = event.oldSize().width()

    def externalResizeDone(self):
        """
        This only happens if the MDI sub window is not maximized and it gets resized; that does
        not normally call the resizeEvent for the alignment window for some reason.
        """
        #print("Done Resizing")
        self.seqArrangeColor()



    def seqInit(self):
        """
        Sequences are stored as triple layer arrays:
        First layer is SEQUENCE
        Second layer is SEQUENCE POSITION
        Third layer is RESIDUE and COLOR
        """
        for seq in self._seqs.values():
            if len(seq) > self.maxlen:
                self.maxlen = len(seq)
        for name, seq in self._seqs.items():
            self.splitNames.append(name)
            local = []
            for i in range(self.maxlen):
                try:
                    char = seq[i]
                    color = self.theme[char]
                    local.append([char, color])
                except IndexError:
                    local.append([" ", None])
                except KeyError:
                    char = seq[i]
                    local.append([char, None])
            self.splitSeqs.append(local)

    def seqArrangeNoColor(self):
        #print("NO COLOR")
        nseqs = len(self._seqs.keys())  # Calculate number of sequences
        charpx = self.fmF.averageCharWidth()
        width = self.alignPane.size().width() - 30
        char_count = int(width / charpx - 20 / charpx)
        if self.alignPane.verticalScrollBar().isVisible():
            char_count = int(width / charpx - 20 / charpx - \
                             self.alignPane.verticalScrollBar().size().width() / charpx)

        lines = int(self.maxlen / char_count)
        self.alignPane.clear()
        self.alignPane.setTextBackgroundColor(Qt.white)
        nline = 0
        for line in range(lines):
            start = nline * char_count
            end = nline * char_count + char_count
            for n in range(nseqs):
                self.alignPane.append("".join([x[0] for x in self.splitSeqs[n][start:end]]))
            self.alignPane.append("")
            self.alignPane.moveCursor(QTextCursor.Start)
            nline += 1

    def seqArrangeColor(self):
        #print("COLOR")
        nseqs = len(self._seqs.keys())  # Calculate number of sequences
        charpx = self.fmF.averageCharWidth()
        width = self.alignPane.size().width() - 30
        char_count = int(width / charpx - 20 / charpx)
        if self.alignPane.verticalScrollBar().isVisible():
            char_count = int(width / charpx - 20 / charpx - \
                             self.alignPane.verticalScrollBar().size().width() / charpx)

        lines = int(self.maxlen/char_count)
        self.alignPane.clear()
        nline = 0
        for line in range(lines):
            start = nline * char_count
            end = nline * char_count + char_count
            for n in range(nseqs):
                for i in range(start, end):
                    self.alignPane.moveCursor(QTextCursor.End)
                    self.alignPane.setTextBackgroundColor(Qt.white)
                    if self.splitSeqs[n][i][1]:
                        self.alignPane.setTextBackgroundColor(self.splitSeqs[n][i][1])
                    self.alignPane.insertPlainText(self.splitSeqs[n][i][0])
                    self.alignPane.setTextBackgroundColor(Qt.white)
                self.alignPane.insertPlainText("\n")
            nline += 1
            self.alignPane.insertPlainText("\n")
        self.alignPane.moveCursor(QTextCursor.Start)

    def seqs(self):
        return self._seqs

    def setSeqs(self, seqs):
        self._seqs = seqs
        self.seqArrangeColor()

    def updateName(self, old, new):
        seq = self._seqs[old]
        self._seqs[new]=seq
        self._seqs.pop(old)
        self.seqArrangeColor()


class QuitDialog(QDialog, quit_ui.Ui_closeConfirm):
    def __init__(self, parent):
        super(self.__class__, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Leaving so soon?")
        self.buttonBox.button(QDialogButtonBox.Ok).clicked.connect(self.ok)
        self.buttonBox.button(QDialogButtonBox.Discard).clicked.connect(self.discard)

    def ok(self):
        self.done(1)

    def discard(self):
        self.done(2)
