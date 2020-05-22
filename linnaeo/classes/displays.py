from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor, QFontDatabase, QFont, QFontMetricsF, QTextCursor, QCursor
from PyQt5.QtWidgets import QWidget, QApplication, QDialog, QDialogButtonBox

from linnaeo.classes import widgets
from linnaeo.resources import linnaeo_rc
from linnaeo.ui import alignment_ui, quit_ui


class AlignSubWindow(QWidget, alignment_ui.Ui_aliWindow):
    """
    Alignment SubWindow UI. Takes in a dictionary of sequences that have been aligned and arranges them.
    # TODO: This does not maintain the alignment order. Shant be helped?...
    """
    resized = pyqtSignal()
    toolTipReq = pyqtSignal()


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
        #self.aliDoc = widgets.AlignDoc()
        #self.alignPane.setDocument(self.aliDoc)
        self._seqs = seqs
        self.resized.connect(self.externalResizeDone)
        self.alignPane.verticalScrollBar().valueChanged.connect(self.namePane.verticalScrollBar().setValue)
        self.toolTipReq.connect(self.getSeqTT)
        self.theme = None
        self.splitNames = []
        self.splitSeqs = []
        self.setMouseTracking(True)

        # FANCY FONTWORK
        # This maintains the font within the application.
        fid = QFontDatabase.addApplicationFont(":/fonts/LiberationMono.ttf")
        family = QFontDatabase.applicationFontFamilies(fid)[0]
        font = QFont(family, 10)
        self.fmF = QFontMetricsF(font)  # FontMetrics Float... because default FontMetrics gives Int
        self.alignPane.document().setDefaultFont(font)
        self.alignPane.setStyleSheet("QTextEdit {padding-left:20px; padding-right:0px; background-color: \
                                     rgb(255,255,255)}")
        self.namePane.setStyleSheet("QTextEdit {padding-top:1px;}")
        self.alignPane.setCursorWidth(0)

        self.refseq = None
        self.maxlen = 0
        self.maxname = 0

        # options to do
        # TODO: Implement these
        self.showRuler = True
        self.showColors = True  # Partially implemented
        self.relColors = False

        if self.showColors:
            self.theme = self.defTheme

        self.seqInit()

    def eventFilter(self, obj, event):
        if event.type() == 5:
            print("MOUSE MOVE")
            return True
        else:
            return super().eventFilter(obj, event)


    def mouseMoveEvent(self, event):
        #print(QCursor.pos())
        return super().mouseMoveEvent(event)

    def event(self, event):
        #print(event, event.type())
        #if event.type() == 110:
        #    self.toolTipReq.emit()
        return super().event(event)

    def getSeqTT(self):
        print("TOOLTIP")
        pos = QCursor.pos()
        self.alignPane.setToolTip("%s, %s" % (str(pos.x()), str(pos.y())))
        self.alignPane.toolTip()


    def setTheme(self, theme):
        self.theme = theme

    def toggleRulers(self):
        self.showRuler = not self.showRuler
        self.seqArrange()

    def toggleColors(self):
        self. showColors = not self.showColors
        self.seqArrange()

    def mousePressEvent(self, event):
        print(self.alignPane.textCursor().positionInBlock())
        print(self.alignPane.textCursor().position())
        return super().mousePressEvent(event)

    def resizeEvent(self, event):
        """
        This gets called anytime the window is in the process of being redrawn. If the MDI Subwindow is maximized,
        it calls a resizeEvent upon release too.
        """
        if self.userIsResizing:
            self.seqArrange(color=False)
        elif not self.userIsResizing:
            self.resized.emit()
        super(AlignSubWindow, self).resizeEvent(event)
        # self.oldwidth = event.oldSize().width()

    def externalResizeDone(self):
        """
        This only happens if the MDI sub window is not maximized and it gets resized; that does
        not normally call the resizeEvent for the alignment window for some reason.
        """
        self.seqArrange()

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
            if len(name) > self.maxname:
                self.maxname = len(name)
            local = []
            count = 0
            for i in range(self.maxlen):
                try:
                    char = seq[i]
                    color = self.theme[char]
                    local.append([char, color, count])
                    if char != "-":
                        count += 1
                except IndexError:
                    local.append([" ", None])
                except KeyError:
                    char = seq[i]
                    local.append([char, None])
            self.splitSeqs.append(local)

    def seqArrange(self, color=True):
        if not self.showColors:
            color = False
        nseqs = len(self._seqs.keys())  # Calculate number of sequences
        charpx = self.fmF.averageCharWidth()
        width = self.alignPane.size().width() - 30
        char_count = int(width / charpx - 20 / charpx)
        if self.alignPane.verticalScrollBar().isVisible():
            char_count = int(width / charpx - 20 / charpx - \
                             self.alignPane.verticalScrollBar().size().width() / charpx)
        lines = int(self.maxlen/char_count)+1
        self.alignPane.clear()
        self.namePane.clear()
        self.namePane.setMinimumWidth((self.maxname * charpx) + 5)
        nline = 0
        for line in range(lines):
            self.alignPane.moveCursor(QTextCursor.End)
            start = nline * char_count
            end = nline * char_count + char_count
            gap = 0
            if line == lines-1:
                oldend = end
                end = start+len(self.splitSeqs[0][start:])
                #gap = oldend - end
            if self.showRuler:
                self.namePane.insertPlainText("\n")
                ruler = str(start+1)+" "*(char_count-gap-len(str(start+1))-len(str(end)))+str(end)
                self.alignPane.insertPlainText(ruler)
                self.alignPane.insertPlainText("\n")
            for n in range(nseqs):
                self.namePane.setAlignment(Qt.AlignRight)
                self.namePane.insertPlainText(self.splitNames[n])
                self.alignPane.moveCursor(QTextCursor.End)
                if not color:
                    self.alignPane.insertPlainText("".join([x[0] for x in self.splitSeqs[n][start:end]]))
                    self.alignPane.insertPlainText("\n")
                elif color:
                    for i in range(start, end):

                        self.alignPane.moveCursor(QTextCursor.End)
                        self.alignPane.setTextBackgroundColor(Qt.white)
                        if self.splitSeqs[n][i][1]:
                            self.alignPane.setTextBackgroundColor(self.splitSeqs[n][i][1])
                        self.alignPane.insertPlainText(self.splitSeqs[n][i][0])
                        self.alignPane.setTextBackgroundColor(Qt.white)
                    self.alignPane.insertPlainText("\n")
                self.namePane.insertPlainText("\n")
            nline += 1
            self.alignPane.insertPlainText("\n")
            self.namePane.insertPlainText("\n")
        self.alignPane.moveCursor(QTextCursor.Start)
        self.namePane.moveCursor(QTextCursor.Start)


    def seqs(self):
        return self._seqs

    def setSeqs(self, seqs):
        self._seqs = seqs
        self.seqArrange()

    def updateName(self, old, new):
        seq = self._seqs[old]
        self._seqs[new]=seq
        self._seqs.pop(old)
        self.seqArrange()


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
