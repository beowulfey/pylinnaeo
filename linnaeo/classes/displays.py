import logging
from math import floor

from PyQt5.QtCore import pyqtSignal, Qt, QUrl
from PyQt5.QtGui import QFontMetricsF, QColor
# from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import QWidget, QDialog, QDialogButtonBox, QPushButton, QMainWindow, QTextEdit, QFrame, QSizePolicy

from linnaeo.classes import widgets, utilities, themes
from linnaeo.ui import alignment_ui, quit_ui, about_ui, ali_settings_ui, comments_ui


class NGLviewer(QMainWindow):
    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        webview = QWebEngineView(self)
        webview.load(QUrl('http://nglviewer.org/ngl/'))
        self.setCentralWidget(webview)


class CommentsPane(QWidget, comments_ui.Ui_Form):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)


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
        self.rulerPane = QTextEdit()
        self.commentPane = CommentsPane()
        self.commentButton = QPushButton("Save")

        # Init functional variables
        self.fmF = None
        self._seqs = seqs
        self.dssps = {}
        self.splitNames = []
        self.splitSeqs = []
        self.userIsResizing = False
        self.refseq = None
        self.last = None
        self.maxlen = 0
        self.maxname = 0
        self.lines = 0
        self.comments = {}
        self.showRuler = False
        self.showColors = False
        self.consvColors = False
        self.showDSSP = True

        # Draw the window
        self.setupUi(self)
        self.setupCustomUi()

        # Connect all slots and start
        self.rulerPane.verticalScrollBar().valueChanged.connect(self.alignPane.verticalScrollBar().setValue)
        self.alignPane.verticalScrollBar().valueChanged.connect(self.rulerPane.verticalScrollBar().setValue)
        self.alignPane.verticalScrollBar().valueChanged.connect(self.namePane.verticalScrollBar().setValue)
        self.namePane.verticalScrollBar().valueChanged.connect(self.alignPane.verticalScrollBar().setValue)
        self.nameChange.connect(self.updateName)
        self.lineChange.connect(self.nameArrange)
        self.alignPane.commentAdded.connect(self.showCommentWindow)

        # Initialize settings
        self.theme = self.lookupTheme('Default')
        self.params = {}
        self.setParams(params)

        self.done = True
        self.seqInit()

    def setupCustomUi(self):
        #self.rulerPane.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        sizePolicy = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.rulerPane.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding))
        self.namePane.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding))
        self.alignPane.setSizePolicy(QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding))
        self.alignPane.setMinimumWidth(100)
        self.rulerPane.setMaximumWidth(50)
        self.rulerPane.setMinimumWidth(5)
        self.alignPane.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.alignPane.setFrameShape(QFrame.StyledPanel)
        self.alignPane.setFrameShadow(QFrame.Plain)
        self.rulerPane.setFrameShape(QFrame.NoFrame)
        self.namePane.setFrameShape(QFrame.NoFrame)
        self.rulerPane.setCursorWidth(0)
        #self.namePane.setLineWidth(0)
        #self.alignPane.setLineWidth(0)
        #self.rulerPane.setLineWidth(0)
        self.alignPane.setStyleSheet("QTextEdit {padding-left:20px; padding-right:0px; background-color: \
                                             rgb(255,255,255)}")
        self.namePane.setStyleSheet("QTextEdit {padding-top:1px;background-color:\
                                                rgb(238, 238, 239)}")
        self.rulerPane.setStyleSheet("QTextEdit {padding-top:1px;padding-left:0px; padding-right:0px; background-color:\
                                                rgb(238, 238, 239)}")
        self.gridLayout_2.addWidget(self.alignPane, 0, 1)
        self.gridLayout_2.addWidget(self.rulerPane, 0, 2)

    def seqInit(self):
        """
        Sequences are stored as double layer arrays of each character.
        First layer is SEQUENCE
        Second layer is SEQUENCE POSITION
        Third layer is RESIDUE and COLOR
        so self.seq = [SEQ = [Position = [HTML Char/Color, ResID, Dssp], ... ], ... ]
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
                    dssp = None
                    color = self.theme[char]
                    if self.theme == themes.Comments().theme:
                        if i in self.comments.keys():
                            print(self.comments[i])
                            color = QColor(Qt.yellow)
                    tcolor = '#FFFFFF' if color.getHsl()[2] / 255 * 100 <= 50 else '#000000'
                    char = '<span style=\"background-color: %s; color: %s\">%s</span>' % (
                        color.name(), tcolor, char)
                    if self.dssps:
                        index = list(self._seqs.values()).index(seq)
                        try:
                            print("adding dssp")
                            dssp = self.dssps[index][tcount]
                        except KeyError:
                            dssp = '-'
                    else:
                        dssp = None
                else:
                    char = '<span style=\"background-color:#FFFFFF;\">' + seq[i] + "</span>"
                    tcount = 0
                    dssp = "-"
                local.append([char, tcount, dssp])
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

    def seqArrange(self, color=True, rulers=True, dssp=True):
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
            if not self.showDSSP:
                dssp = False
            charpx = self.fmF.averageCharWidth()
            self.rulerPane.size().setWidth(4 * charpx + 3)
            width = self.alignPane.size().width()
            print(width)
            char_count = int(width / charpx - 40 / charpx)
            if self.rulerPane.verticalScrollBar().isVisible():
                self.rulerPane.resize(int(4 * charpx + 3 + (self.rulerPane.verticalScrollBar().size().width())),
                                      self.rulerPane.size().height())
                sb = self.rulerPane.verticalScrollBar()
                self.last = self.rulerPane.verticalScrollBar().sliderPosition() / (
                            self.rulerPane.verticalScrollBar().maximum() -
                            self.rulerPane.verticalScrollBar().minimum())
            lines = int(self.maxlen / char_count) + 1
            if lines != self.lines:
                self.lineChange.emit(lines)
                self.lines = lines
            self.alignPane.lines = lines
            self.alignPane.setChars(char_count)
            self.alignPane.names = self.splitNames
            self.alignPane.clear()
            fancy = False if self.userIsResizing else True
            worker = utilities.SeqThread(self.splitSeqs, char_count, lines, rulers, color, dssp, fancy=fancy, parent=self, )
            worker.start()
            worker.wait()
            style = "<style>pre{font-family:%s; font-size:%spt;}</style>" % (
                self.font().family(), self.font().pointSize())
            self.alignPane.setHtml(style + worker.html)
            # RULER CALCULATION --> SIDE PANEL.
            self.rulerPane.clear()
            rulerHtml = ["<pre style=\"font-family:%s; font-size:%spt; text-align: left;\">" %
                         (self.font().family(),self.font().pointSize())]
            for x in range(self.lines):
                if self.showRuler and self.showDSSP:
                    exline = "\n\n\n"
                elif self.showRuler or self.showDSSP:
                    exline = "\n\n"
                else:
                    exline = "\n"
                #exline = "\n\n" if self.showRuler else "\n"
                rulerHtml.append(exline)
                for i in range(len(self.splitSeqs)):
                    try:
                        label = str(self.splitSeqs[i][x * char_count + char_count - 1][1])
                        if label == "0":
                            for y in range(char_count):
                                label = str(self.splitSeqs[i][x * char_count + char_count - 1 - y][1])
                                if label != "0":
                                    break
                    except IndexError:
                        label = ""
                    rulerHtml.append(label + "\n")
            rulerHtml.append('</pre>')
            self.rulerPane.setHtml(style + "".join(rulerHtml))
            if self.rulerPane.verticalScrollBar().isVisible():
                if self.last:
                    self.rulerPane.verticalScrollBar().setSliderPosition(int(round(((self.rulerPane.verticalScrollBar().maximum() -
                                                                            self.rulerPane.verticalScrollBar().minimum()) *
                                                                           self.last))))
        except ZeroDivisionError:
            self.alignLogger.info("Font returned zero char width. Please choose a different font")

    # UTILITY FUNCTIONS
    def setTheme(self, theme):
        self.theme = self.lookupTheme(theme)
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
            # print("IGNORING")
            font.setPointSize(self.font().pointSize())
        # print(font.pointSize())
        super().setFont(font)
        # print("FINAL", self.font().pointSize())
        self.fmF = QFontMetricsF(self.font())
        if self.done:
            # print("REDRAWING")
            self.seqInit()
            self.seqArrange()
            self.nameArrange(self.lines)

    def setFontSize(self, size):
        font = self.font()
        font.setPointSize(size)
        # print("FONT",font.pointSize(), self.font().pointSize())
        self.setFont(font)

    def setParams(self, params):
        # print("UPDATING VALUES")
        self.params = params
        # print(self.params)
        self.showRuler = self.params['ruler']
        self.showColors = self.params['colors']
        self.consvColors = self.params['byconsv']
        if self.font().pointSize() != self.params['fontsize']:
            # print("Changing font size")
            self.setFontSize(self.params['fontsize'])
        if self.font() != self.params['font']:
            self.setFont(self.params['font'])
        newtheme = self.lookupTheme(self.params['theme'])
        if self.theme != newtheme:
            self.theme = newtheme
            if self.done:
                self.seqInit()
                self.seqArrange()

    def lookupTheme(self, theme):
        """ Converts the stored theme name into a class """
        match = themes.PaleByType().theme
        if theme == 'Default':
            match = themes.PaleByType().theme
        elif theme == 'Bold':
            match = themes.Bold().theme
        elif theme == 'Monochrome':
            match = themes.Mono().theme
        elif theme == 'ColorSafe':
            match = themes.ColorSafe().theme
        elif theme == 'Rainbow':
            match = themes.Rainbow().theme
        elif theme == 'Grayscale':
            match = themes.Grayscale().theme
        elif theme == 'Annotations':
            match = themes.Comments().theme
        return match

    def showCommentWindow(self, target):
        # TODO: shows for ALL rows!
        name = self.splitNames[target[0]]
        resi = self.splitSeqs[target[0]][target[2]]
        self.comments[target[2]] = "COMMENT"
        self.seqInit()
        self.seqArrange()
        self.commentPane.lineEdit.setText(str(name) + " " + str(resi))
        # self.gridLayout.addWidget(self.commentButton,1,0)
        self.gridLayout.addWidget(self.commentPane, 1, 1)

    def addStructure(self, dssp, seq):
        """ Adds DSSP data to the SplitSeqs array. """
        seqs = list(self._seqs.values())
        if str(seq.seq) in seqs:
            index = seqs.index(str(seq.seq))
            self.dssps[index]=dssp
            for res in self.splitSeqs[index]:
                if not res[2]:
                    try:
                        res[2] = dssp[res[1]]
                    except KeyError:
                        res[2] = "-"
                else:
                    print("Weird, got a duplicate at %s " % res[1])
            print(self.splitSeqs[index])




class OptionsPane(QWidget, ali_settings_ui.Ui_Form):
    # updateParam = pyqtSignal(dict)

    def __init__(self, parent):
        super().__init__(parent)
        self.setupUi(self)
        self.setFixedWidth(140)
        self.params = {}
        self.themeIndices = {}
        self.initPane()

    def initPane(self):
        for index in range(0, self.comboTheme.model().rowCount()):
            self.themeIndices[self.comboTheme.model().itemData(self.comboTheme.model().index(index, 0))[0]] = index

        self.checkRuler.toggled.connect(self.rulerToggle)
        self.checkColors.toggled.connect(self.colorToggle)
        self.comboTheme.currentIndexChanged.connect(self.changeTheme)
        self.comboFont.currentFontChanged.connect(self.changeFont)
        self.spinFontSize.valueChanged.connect(self.changeFontSize)

    def setParams(self, params):
        """ These are set by the preferences pane --> default for every new window """
        self.params = params.copy()
        #print(params["font"].family(),
         #     self.params["font"].family())  # 'rulers', 'colors', 'fontsize', 'theme', 'font', 'byconsv'
        self.checkRuler.setChecked(self.params['ruler'])
        self.checkColors.setChecked(self.params['colors'])
        self.checkConsv.setChecked(self.params['byconsv'])
        self.comboTheme.setCurrentIndex(self.themeIndices[self.params['theme']])
        self.spinFontSize.setValue(self.params['fontsize'])
        self.comboFont.setCurrentFont(self.params['font'])
        #print(self.params['font'].family())

    def rulerToggle(self):
        self.params['ruler'] = self.checkRuler.isChecked()

    def colorToggle(self):
        self.params['colors'] = self.checkColors.isChecked()

    def changeTheme(self):
        self.params['theme'] = self.comboTheme.currentText()

    def changeFont(self):
        self.params['font'] = self.comboFont.currentFont()

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
