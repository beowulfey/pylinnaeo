from PyQt5.QtCore import pyqtSignal
from PyQt5.QtGui import QFontDatabase, QFont, QFontMetricsF
from PyQt5.QtWidgets import QWidget, QDialog, QDialogButtonBox

from linnaeo.classes import widgets, utilities, themes
from linnaeo.ui import alignment_ui, quit_ui, about_ui

from linnaeo.resources import linnaeo_rc
class AlignSubWindow(QWidget, alignment_ui.Ui_aliWindow):
    """
    Alignment SubWindow UI. Takes in a dictionary of sequences that have been aligned and arranges them.
    Sequences in the alignment are threaded and prepared as HTML; the color is stored in the array in order to
    reduce the load on display. However, just drawing the display is computationally expensive (as is generating
    the ruler for it), so both are turned off during resizing.
    # TODO: Currently, the shut off does not happen on MacOS because Qt on macs doesn't alert on window presses. Damn.
    # TODO: This does not maintain the alignment order -- trouble with the ClustalO API. Shant be helped?...
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
        self.fmF = QFontMetricsF(self.font)  # FontMetrics Float... because default FontMetrics gives Int, which is bad.

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


        # Draw the window
        self.setupUi(self)
        self.setupCustomUi()

        # Connect all slots and start
        self.resized.connect(self.externalResizeDone)
        self.alignPane.verticalScrollBar().valueChanged.connect(self.namePane.verticalScrollBar().setValue)
        self.namePane.verticalScrollBar().valueChanged.connect(self.alignPane.verticalScrollBar().setValue)
        self.nameChange.connect(self.updateName)
        self.lineChange.connect(self.nameArrange)

        self.seqInit()


    def setupCustomUi(self):
        self.horizontalLayout.addWidget(self.alignPane)
        self.alignPane.document().setDefaultFont(self.font)
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
            if len(seq) < self.maxlen:
                for n in range(self.maxlen - len(seq)):
                    seq.append(" ")
            local = []
            count = 0
            for i in range(self.maxlen):
                color = None
                char = seq[i]
                if char not in ["-", " "]:
                    count += 1
                    tcount = count
                    color = self.theme[char]
                    char = '<span style=\"background-color:'+color.name()+';\">'+char+"</span>"
                else:
                    char = '<span style=\"background-color:#FFFFFF;\">' + seq[i] + "</span>"
                    tcount = 0
                local.append([char, tcount])
            self.splitSeqs.append(local)
        #self.seqArrange()
        self.alignPane.seqs = self.splitSeqs

    def nameArrange(self, lines):
        """ Generates the name panel; only fires if the number of lines changes to avoid needless computation"""
        self.namePane.clear()
        self.namePane.setMinimumWidth((self.maxname * self.fmF.averageCharWidth()) + 5)
        names = ["<pre style=\"text-align: right;\">\n"]
        for line in range(lines):
            if self.showRuler:
                names.append("\n")
            for i in range(len(self.splitNames)):
                names.append(self.splitNames[i] + "\n")
            names.append("\n")
        names.append("</pre>")
        self.namePane.setHtml("".join(names))

    def seqArrange(self, color=True, rulers=True):
        """
        The bread and butter. This fires upon creation and any resizing event. Computing the resize is done
        in a separate thread to help with smoothness; showing color and rulers is still very slow though. Resize events
        call this function with color off, and the ruler is turned off automatically.
        """
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
        self.alignPane.setChars(char_count)
        self.alignPane.names = self.splitNames
        if lines != self.lines:
            self.lineChange.emit(lines)
            self.lines = lines
        self.alignPane.lines = lines
        self.alignPane.clear()
        fancy = False if self.userIsResizing else True
        worker = utilities.SeqThread(self.splitSeqs, char_count, lines, rulers, color, fancy)
        worker.start()
        worker.wait()
        self.alignPane.setText(worker.html)

    # UTILITY FUNCTIONS
    def setTheme(self, theme):
        self.theme = theme

    def toggleRulers(self):
        self.showRuler = not self.showRuler
        self.nameArrange(self.lines)
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
        self.alignPane.setFont(self.font)
        self.namePane.setFont(self.font)
        self.nameArrange(self.lines)
        self.seqArrange()

    def decreaseFont(self):
        size = self.font.pointSizeF() - 1
        self.font.setPointSizeF(size)
        self.fmF = QFontMetricsF(self.font)
        self.namePane.setFont(self.font)
        self.nameArrange(self.lines)
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
