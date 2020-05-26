from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor, QFontDatabase, QFont, QFontMetricsF, QTextCursor, QCursor
from PyQt5.QtWidgets import QWidget, QApplication, QDialog, QDialogButtonBox

from linnaeo.classes import widgets, utilities, themes
from linnaeo.resources import linnaeo_rc
from linnaeo.ui import alignment_ui, quit_ui, about_ui


class AlignSubWindow(QWidget, alignment_ui.Ui_aliWindow):
    """
    Alignment SubWindow UI. Takes in a dictionary of sequences that have been aligned and arranges them.
    # TODO: This does not maintain the alignment order. Shant be helped?...
    """
    resized = pyqtSignal()
    nameChange = pyqtSignal((str, str))
    lineChange = pyqtSignal(int)

    def __init__(self, seqs):
        super(self.__class__, self).__init__()

        # Construct the window
        self.alignPane = widgets.AlignPane(self)
        self.family = (QFontDatabase.applicationFontFamilies(QFontDatabase.addApplicationFont(
            ':/fonts/LiberationMono.ttf'))[0])
        self.font = QFont(self.family, 10)
        self.fmF = QFontMetricsF(self.font) # FontMetrics Float... because default FontMetrics gives Int

        # Initialize settings
        self.theme = themes.PaleTheme().theme
        self.showRuler = True
        self.showColors = True
        self.relColors = False

        # Init functional variables
        self._seqs = seqs
        self.splitNames = []
        self.splitSeqs = []
        self.userIsResizing = False
        self.refseq = None
        self.lastpos = None
        self.maxlen = 0
        self.maxname = 0
        self.lines = 0

        self.seqInit()

        # Draw the window
        self.setupUi(self)
        self.setupCustomUi()

        # Connect all slots and start
        self.resized.connect(self.externalResizeDone)
        self.alignPane.verticalScrollBar().valueChanged.connect(self.namePane.verticalScrollBar().setValue)
        self.nameChange.connect(self.updateName)
        self.lineChange.connect(self.nameArrange)


    def setupCustomUi(self):
        self.horizontalLayout.addWidget(self.alignPane)
        self.alignPane.setFont(self.font)
        self.alignPane.setStyleSheet("QTextEdit {padding-left:20px; padding-right:0px; background-color: \
                                             rgb(255,255,255)}")
        self.namePane.setStyleSheet("QTextEdit {padding-top:1px;}")

    def seqInit(self):
        """
        Sequences are stored as double layer arrays of each character.
        First layer is SEQUENCE
        Second layer is SEQUENCE POSITION
        Third layer is RESIDUE and COLOR
        so self.seq = [SEQ = [Position = [HTML Char/Color, ResID], ... ], ... ]
        """
        self.splitSeqs = []
        self.splitNames = []
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
                    char = '<span style=\"background-color:'+color.name()+';\">'+char+"</span>"
                    local.append([char, count])
                    if char != "-":
                        count += 1
                except IndexError:
                    print("NoIndex")
                    local.append([" ", None])
                except KeyError:
                    print("noKey")
                    char = '<span style=\"background-color:#FFFFFF;\">'+seq[i]+"</span>"
                    local.append([char, None])
            self.splitSeqs.append(local)

    def nameArrange(self, lines):
        self.namePane.clear()
        self.namePane.setMinimumWidth((self.maxname * self.fmF.averageCharWidth()) + 5)
        names = ["<pre style=\"text-align: right;\">"]
        if self.showRuler:
            names.append("\n")
        for line in range(lines):
            for i in range(len(self.splitNames)):
                if self.showRuler:
                    names.append("\n")
                names.append(self.splitNames[i] + "\n")
            names.append("\n")
        names.append("</pre>")
        self.namePane.setHtml("".join(names))

    def seqArrange(self, color=True, rulers=True):
        if not self.showColors:
            color = False
        if not self.showRuler:
            rulers = False
        # Calculate font and window metrics
        charpx = self.fmF.averageCharWidth()
        width = self.alignPane.size().width() - 30
        char_count = int(width / charpx - 20 / charpx)
        if self.alignPane.verticalScrollBar().isVisible():
            char_count = int(width / charpx - 20 / charpx - \
                             self.alignPane.verticalScrollBar().size().width() / charpx)
        lines = int(self.maxlen / char_count) + 1
        if lines != self.lines:
            self.lineChange.emit(lines)
        self.alignPane.clear()
        if self.userIsResizing:
            self.alignPane.setText(utilities.redrawBasic(self.splitSeqs, char_count, lines, rulers))
        else:
            self.alignPane.setHtml(utilities.redrawFancy(self.splitSeqs, char_count, lines, rulers, color))

    def seqArrangeBk(self, color=True, rulers=True):
        if not self.showColors:
            color = False
        if not self.showRuler:
            rulers = False
        # Calculate font and window metrics
        nseqs = len(self._seqs.keys())
        charpx = self.fmF.averageCharWidth()
        width = self.alignPane.size().width() - 30
        char_count = int(width / charpx - 20 / charpx)
        if self.alignPane.verticalScrollBar().isVisible():
            char_count = int(width / charpx - 20 / charpx - \
                             self.alignPane.verticalScrollBar().size().width() / charpx)
        lines = int(self.maxlen / char_count) + 1
        self.namePane.setMinimumWidth((self.maxname * charpx) + 5)
        nline = 0
        # if not self.userIsResizing:
        # TODO: ISOLATE THIS OUT TO DO IN THREADS
        self.alignPane.clear()
        self.namePane.clear()
        for line in range(lines):
            self.alignPane.moveCursor(QTextCursor.End)
            start = nline * char_count
            end = nline * char_count + char_count
            gap = 0
            if line == lines - 1:
                oldend = end
                end = start + len(self.splitSeqs[0][start:])
                gap = oldend - end
            if self.showRuler and rulers:  # TODO REMOVE AND NOT
                self.namePane.insertPlainText("\n")
                ruler = utilities.buildRuler(char_count, gap, start, end)
                self.alignPane.insertPlainText(ruler)
                self.alignPane.insertPlainText("\n")
            elif self.showRuler and not rulers:  # TODO MOVE THIS TO BELOW
                self.namePane.insertPlainText("\n")
                self.alignPane.insertPlainText(" ")
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
        if self.lastpos:
            self.alignPane.moveCursor(self.lastpos)
            self.namePane.moveCursor(self.lastpos)
        else:
            self.alignPane.moveCursor(QTextCursor.Start)
            self.namePane.moveCursor(QTextCursor.Start)

    # elif self.userIsResizing:  # TODO DO IN SEPARATE THREAD
    #    print("RESIZING")
    #    self.alignPane.clear()
    #    self.namePane.clear()



    # UTILITY FUNCTIONS
    def setTheme(self, theme):
        self.theme = theme

    def toggleRulers(self):
        self.showRuler = not self.showRuler
        self.seqArrange()

    def toggleColors(self):
        self.showColors = not self.showColors
        self.seqArrange()

    def resizeEvent(self, event):
        """
        This gets called anytime the window is in the process of being redrawn. If the MDI Subwindow is maximized,
        it calls a resizeEvent upon release too.
        """
        if self.userIsResizing:
            self.seqArrange(color=False) #rulers=False)
        elif not self.userIsResizing:
            self.resized.emit()
        super(AlignSubWindow, self).resizeEvent(event)

    def externalResizeDone(self):
        """
        This only happens if the MDI sub window is not maximized and it gets resized; that does
        not normally call the resizeEvent for the alignment window for some reason.
        """
        self.seqArrange()

    def seqs(self):
        return self._seqs

    def setSeqs(self, seqs):
        self._seqs = seqs
        self.seqArrange()

    def updateName(self, old, new):
        seq = self._seqs[old]
        self._seqs[new] = seq
        self._seqs.pop(old)
        self.seqInit()
        self.seqArrange()

    def increaseFont(self):
        size = self.font.pointSizeF() + 1
        self.font.setPointSizeF(size)
        self.fmF = QFontMetricsF(self.font)
        print(self.fmF.averageCharWidth())
        self.alignPane.setFont(self.font)
        self.seqArrange()

    def decreaseFont(self):
        size = self.font.pointSizeF() - 1
        self.font.setPointSizeF(size)
        self.fmF = QFontMetricsF(self.font)
        print(self.fmF.averageCharWidth())
        self.alignPane.setFont(self.font)
        self.seqArrange()

    def getFontSize(self):
        return self.font.pointSize()


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


class AboutDialog(QDialog, about_ui.Ui_Dialog):
    def __init__(self, parent):
        super(self.__class__, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Thanks for reading this!")
        self.buttonBox.button(QDialogButtonBox.Ok).clicked.connect(self.ok)

    def ok(self):
        self.done(1)
