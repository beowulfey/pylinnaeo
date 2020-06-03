import logging

from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QFontDatabase, QFont, QFontMetricsF, QStandardItem, QColor
from PyQt5.QtWidgets import QWidget, QDialog, QDialogButtonBox, qApp

from linnaeo.classes import widgets, utilities, themes
from linnaeo.ui import alignment_ui, quit_ui, about_ui, ali_settings_ui
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

    def __init__(self, seqs, params):
        super(self.__class__, self).__init__()
        self.done = False
        # Construct the window
        self.alignLogger = logging.getLogger("AlignWindow")
        self.alignPane = widgets.AlignPane(self)

        # Init functional variables
        self._seqs = seqs
        self.splitNames = []
        self.splitSeqs = []
        self.userIsResizing = False
        self.refseq = None
        self.last = None
        self.maxlen = 0
        self.maxname = 0
        self.lines = 0

        # Draw the window
        self.setupUi(self)
        self.setupCustomUi()

        # Connect all slots and start
        #self.resized.connect(self.externalResizeDone)
        self.alignPane.verticalScrollBar().valueChanged.connect(self.namePane.verticalScrollBar().setValue)
        self.namePane.verticalScrollBar().valueChanged.connect(self.alignPane.verticalScrollBar().setValue)
        self.nameChange.connect(self.updateName)
        self.lineChange.connect(self.nameArrange)

        # Initialize settings
        self.theme = self.convertTheme('Default')
        self.params = {}
        self.setParams(params)

        #self.fmF = QFontMetricsF(self.font())  # FontMetrics Float...
        #self.setFont(QFont(self.params['font']))
        #print("AFTER LOAD",self.font().pointSize())
        self.done = True
        self.seqInit()


    def setupCustomUi(self):
        self.horizontalLayout.addWidget(self.alignPane)
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
        self.maxname = 0
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
                char = seq[i]
                if char not in ["-", " "]:
                    count += 1
                    tcount = count
                    color = self.theme[char]
                    #print(color.name(),color.getHsl()[2]/255*100,(100-(color.getHsl()[2]/255*100))*5)
                    tcolor = '#FFFFFF' if color.getHsl()[2]/255*100 <= 50 else '#000000'
                    char = '<span style=\"background-color: %s; color: %s\">%s</span>' % (color.name(), tcolor, char)
                else:
                    char = '<span style=\"background-color:#FFFFFF;\">' + seq[i] + "</span>"
                    tcount = 0
                local.append([char, tcount])
            self.splitSeqs.append(local)
        self.alignPane.seqs = self.splitSeqs

    def nameArrange(self, lines):
        """ Generates the name panel; only fires if the number of lines changes to avoid needless computation"""
        self.namePane.clear()
        self.namePane.setMinimumWidth((self.maxname * self.fmF.averageCharWidth()) + 5)
        names = ["<pre style=\"font-family:%s; font-size:%spt; text-align: right;\">\n" % (self.font().family(),
                                                                                           self.font().pointSize())]
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
        try:
            self.last = None
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
                # This is for saving the scroll position
                self.last = self.alignPane.verticalScrollBar().value()
            lines = int(self.maxlen / char_count) + 1
            if lines != self.lines:
                self.lineChange.emit(lines)
                self.lines = lines
            self.alignPane.lines = lines
            self.alignPane.setChars(char_count)
            self.alignPane.names = self.splitNames
            self.alignPane.clear()
            fancy = False if self.userIsResizing else True
            worker = utilities.SeqThread(self.splitSeqs, char_count, lines, rulers, color, fancy=fancy)
            worker.start()
            worker.wait()
            style = "<style>pre{font-family:%s; font-size:%spt;}</style>" % (self.font().family(), self.font().pointSize())
            self.alignPane.setHtml(style+worker.html)

            if self.last:
                self.alignPane.verticalScrollBar().setValue(self.last)
        except ZeroDivisionError:
            self.alignLogger.info("Font returned zero char width. Please choose a different font")

    # UTILITY FUNCTIONS
    def setTheme(self, theme):
        self.theme = self.convertTheme(theme)
        self.params['theme'] = theme
        self.seqInit()
        self.seqArrange()

    def toggleRuler(self, state):
        self.showRuler = state
        self.params['ruler'] = state
        self.nameArrange(self.lines)
        self.seqArrange()

    def toggleColors(self, state):
        self.showColors = state
        self.params['colors'] = state
        self.seqArrange()

    def seqs(self):
        return self._seqs

    def setSeqs(self, seqs):
        self._seqs = seqs
        self.seqArrange()

    def updateName(self, old, new):
        self.alignLogger.debug("Received name change alert")
        seq = self._seqs[old]
        self._seqs[new] = seq
        self._seqs.pop(old)
        if self.done:
            self.seqInit()
            self.seqArrange()
            self.nameArrange(self.lines)

    def setFont(self, font):
        # Choosing a new font has a built in size, which is annoying

        if font.family() != self.font().family() and font.pointSize() != self.font().pointSize():
            #print("IGNORING")
            font.setPointSize(self.font().pointSize())
        #print(font.pointSize())
        super().setFont(font)
        #print("FINAL", self.font().pointSize())
        self.fmF = QFontMetricsF(self.font())
        if self.done:
            #print("REDRAWING")
            self.seqInit()
            self.seqArrange()
            self.nameArrange(self.lines)

    def setFontSize(self, size):
        font = self.font()
        font.setPointSize(size)
        #print("FONT",font.pointSize(), self.font().pointSize())
        self.setFont(font)

    def setParams(self, params):
        #print("UPDATING VALUES")
        self.params = params
        #print(self.params)
        self.showRuler = self.params['ruler']
        self.showColors = self.params['colors']
        self.consvColors = self.params['byconsv']
        if self.font().pointSize() != self.params['fontsize']:
            #print("Changing font size")
            self.setFontSize(self.params['fontsize'])
        if self.font() != self.params['font']:
            self.setFont(self.params['font'])
        newtheme = self.convertTheme(self.params['theme'])
        if self.theme != newtheme:
            self.theme = newtheme
            self.seqInit()
            self.seqArrange()

    def convertTheme(self, theme):
        """ Converts the stored theme name into a class """
        if theme == 'Default':
            converted = themes.PaleByType().theme
        elif theme == 'Bold':
            converted = themes.Bold().theme
        elif theme == 'Monochrome':
            converted = themes.Mono().theme
        elif theme == 'ColorSafe':
            converted = themes.ColorSafe().theme
        elif theme == 'Rainbow':
            converted = themes.Rainbow().theme
        elif theme == 'Grayscale':
            converted = themes.Grayscale().theme
        return converted


class OptionsPane(QWidget, ali_settings_ui.Ui_Form):
    #updateParam = pyqtSignal(dict)

    def __init__(self, parent):
        super().__init__(parent)
        self.setupUi(self)
        self.setFixedWidth(140)
        self.params = {}
        self.themeIndices = {}
        self.initPane()

    def initPane(self):
        for index in range(0, self.comboTheme.model().rowCount()):
            self.themeIndices[self.comboTheme.model().itemData(self.comboTheme.model().index(index,0))[0]] = index

        self.checkRuler.toggled.connect(self.rulerToggle)
        self.checkColors.toggled.connect(self.colorToggle)
        self.comboTheme.currentIndexChanged.connect(self.changeTheme)
        self.comboFont.currentFontChanged.connect(self.changeFont)
        self.spinFontSize.valueChanged.connect(self.changeFontSize)

    def setParams(self, params):
        """ These are set by the preferences pane --> default for every new window """
        self.params = params.copy()
        # 'rulers', 'colors', 'fontsize', 'theme', 'font', 'byconsv'
        self.checkRuler.setChecked(self.params['ruler'])
        self.checkColors.setChecked(self.params['colors'])
        self.checkConsv.setChecked(self.params['byconsv'])
        self.comboTheme.setCurrentIndex(self.themeIndices[self.params['theme']])
        self.spinFontSize.setValue(self.params['fontsize'])
        self.comboFont.setCurrentFont(self.params['font'])

    def rulerToggle(self):
        self.params['ruler']=self.checkRuler.isChecked()

    def colorToggle(self):
        self.params['colors'] = self.checkColors.isChecked()

    def changeTheme(self):
        self.params['theme']=self.comboTheme.currentText()

    def changeFont(self):
        self.params['font']=self.comboFont.currentFont()

    def changeFontSize(self):
        self.params['fontsize'] = self.spinFontSize.value()

    def showColorDesc(self):
        self.params['colordesc'] = self.checkColorDesc.isChecked()





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
