#!/usr/bin/env python3

# Bioscience components
import logging
import os
import sys
import time
#import random

#import objgraph
import psutil
# PyQt components
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from PyQt5.QtCore import Qt, QThreadPool, QFile, QIODevice, QDataStream, QDir
from PyQt5.QtGui import QStandardItem, QFontDatabase, QFont, QIcon, QTextCursor, QColor, QFontMetrics, QPalette
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QAbstractItemView, qApp, QWidget, QSizePolicy, \
    QFileDialog, QTextEdit, QTextBrowser

#from pympler import tracker, refbrowser  # summary, muppy

import linnaeo
from linnaeo.classes.utilities import OutputWrapper
from linnaeo.resources import linnaeo_rc
from linnaeo.classes import widgets, utilities, methods, models, displays
from linnaeo.classes.displays import QuitDialog, AlignSubWindow
from linnaeo.ui import linnaeo_ui


class Linnaeo(QMainWindow, methods.Slots, methods.Debug, linnaeo_ui.Ui_MainWindow):
    """
    Constructor for the Main Window of the Linnaeo App
    Contains all the user interface functions, as well as underlying code for major features.
    The bulk of the program is located here. Please see classes.methods for additional methods for this class.
    The other major functionality, drawing the alignments themselves, is found under the displays.AlignSubWindow class.
    """

    # sendParams = pyqtSignal(dict)

    def __init__(self, trees=None, data=None):
        super(self.__class__, self).__init__()
        self.geometry().size().setWidth(1200)
        self.start = time.perf_counter()

        # Initialize UI
        self.setAttribute(Qt.WA_QuitOnClose)
        self.setupUi(self)

        if sys.platform in ['darwin']:#, 'linux']:
            # MacOS uses a windowed format for the .apps -- this is nice and portable, but
            # lacks the terminal for information. So I'm integrating the terminal output into
            # a QTextEdit at the bottom of the screen.
            # TODO: Make this a setting to show
            self.console = QTextEdit()
            fmF = QFontMetrics(self.console.document().defaultFont())
            ht = fmF.height()
            self.console.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            self.console.setMaximumHeight(ht*4)
            del fmF, ht
            self.console.setTextInteractionFlags(Qt.TextSelectableByKeyboard | Qt.TextSelectableByMouse)
            self.gridLayout.addWidget(self.console, 1, 0)
            # self.stdout = OutputWrapper(self, True)
            # self.stdout.outputWritten.connect(self.handleOutput)
            self.stderr = OutputWrapper(self, False)
            self.stderr.outputWritten.connect(self.handleOutput)
            logging.basicConfig(force=True, stream=self.stderr)

        # System instants; status bar and threads
        self.memLabel = QLabel()
        self.mainProcess = psutil.Process(os.getpid())
        self.processTimer = utilities.ProcTimerThread(self)

        self.mainLogger = logging.getLogger("Main")
        self.threadpool = QThreadPool()
        self.mainLogger.info("Threading with a maximum of %d threads" % self.threadpool.maxThreadCount())

        # Project instants and inherent variables for logic.
        self.params = {}
        self.localtime = 0
        self.SequenceRole = Qt.UserRole + 2
        self.WindowRole = Qt.UserRole + 3
        self.StructureRole = Qt.UserRole + 4
        self.lastClickedTree = None
        self.lastAlignment = None
        self.windows = None  # Windows stored as { windex : MDISubWindow }
        self.windex = None  # Acts as identifier for tracking alignments (max 2.1 billion)
        self.sequences = None  # Stored as {WindowID:[SeqRecord(s)] }
        self.titles = None  # maintains a list of sequence titles to confirm uniqueness
        self.bioRoot = None
        self.bioModel = None
        self.projectRoot = None
        self.projectModel = None
        self._currentWindow = None
        self.runningDSSP = []

        # MDI Window
        self.mdiArea = widgets.MDIArea(self)
        self.optionsPane = displays.OptionsPane(self)
        self.gridLayout_2.addWidget(self.mdiArea)
        self.gridLayout_2.addWidget(self.optionsPane, 0, 2)
        self.optionsPane.hide()
        self.colorPane = QTextBrowser()
        self.colorPane.setOpenExternalLinks(True)
        self.colorPane.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.colorPane.setTextInteractionFlags(Qt.LinksAccessibleByMouse)#Qt.NoTextInteraction)
        bgcolor = self.palette().color(self.backgroundRole())
        self.colorPane.setStyleSheet("QTextEdit { border:0;background-color:%s;}" % bgcolor.name())
        del bgcolor
        #

        # Tree stuff
        self.bioTree = widgets.TreeView()
        self.projectTree = widgets.TreeView()
        self.mainLogger.debug("Finished initializing, took %f seconds" % float(time.perf_counter() - self.start))

        # Other functions
        self.guiSet(trees, data)
        self.mainLogger.debug("Gui initializing complete, took %f seconds" % float(time.perf_counter() - self.start))
        self.guiFinalize()
        self.mainLogger.debug("Finalized gui, took %f seconds" % float(time.perf_counter() - self.start))
        self.connectSlots()
        #self.mainLogger.debug("Slots connected after %f seconds" % float(time.perf_counter() - self.start))
        #self.mainLogger.debug("Setup took %f seconds" % float(time.perf_counter() - self.start))
        self.setMouseTracking(True)
        del trees, data

    def handleOutput(self, text, stdout):
        color = self.console.textColor()
        newcolor = color.name() if stdout else QColor(Qt.darkRed).name()
        text = '<span style=\"white-space: pre-wrap; color: %s\">%s</span>' % (
            newcolor, text)
        self.console.moveCursor(QTextCursor.End)
        self.console.insertHtml(text)
        self.console.moveCursor(QTextCursor.PreviousCharacter)
        del color, newcolor, text

    def guiSet(self, trees=None, data=None):
        """ Initialize GUI with default parameters. """
        self.lastClickedTree = None
        self.lastAlignment = {}
        self.windows = {}  # Windows stored as { windex : MDISubWindow }
        self.windex = 0  # Acts as identifier for tracking alignments (max 2.1 billion)
        self.sequences = {}  # Stored as {WindowID:[SeqRecord(s)] }
        self.titles = []  # maintains a list of sequence titles to confirm uniqueness
        #self.mainLogger.debug("guiSet took took %f seconds" % float(time.perf_counter() - self.start))

        # Load default options for windows (from parameters file if saved)
        # if PARAMETERS FILE:
        #   params = FROMFILE
        # print(qApp.instance().defFont.family())
        self.default_params = {'ruler': True, 'colors': True, 'fontsize': 10,
                               'theme': 'Default', 'font': qApp.instance().defFont,
                               'byconsv': False, 'tabbed': False,
                               'darkmode': False, 'dssp': False,
                               }
        
        self.params = self.default_params.copy()
        
        if sys.platform in ['win32', 'darwin']:
            self.params['fontsize'] = 12
        # print(self.params['font'].family())
        self.optionsPane.setParams(self.params)

        # This is fired upon loading a saved workspace.
        if trees:
            self.mainLogger.info("Loading saved workspace!")
            self.bioModel, self.projectModel = trees
            self.sequences, self.titles, self.windex = data
            self.windex = int(self.windex)
            self.bioRoot = self.bioModel.invisibleRootItem().child(0)  # TODO: ELIMINATE USE OF ROOT. Use InvsRoot
            self.projectRoot = self.projectModel.invisibleRootItem().child(0)
            self.windows = {}
            self.rebuildTrees()

        else:
            self.bioRoot = QStandardItem("Folder")
            self.bioModel = widgets.ItemModel(self.windows, seqTree=True)
            self.bioRoot.setData("Folder")
            self.bioModel.appendRow(self.bioRoot)
            self.projectRoot = QStandardItem("Folder")
            self.projectModel = widgets.ItemModel(self.windows)
            self.projectModel.appendRow(self.projectRoot)
            self.mainLogger.debug("After Tree Setup")

        self.bioTree.setModel(self.bioModel)
        self.projectTree.setModel(self.projectModel)
        self.bioModel.setHorizontalHeaderLabels(["Sequences"])
        self.projectModel.setHorizontalHeaderLabels(["Alignments"])
        for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            self.bioTree.setExpanded(node.index(), True)
        for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
            self.projectTree.setExpanded(node.index(), True)
        # self.installEventFilter(self)
        del trees, data, node

    def guiFinalize(self):
        # Tree setup
        #self.mainLogger.debug("Up to selectionMode took %f seconds" % float(time.perf_counter() - self.start))
        self.bioTree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        #self.mainLogger.debug("Up to adding bioTree took %f seconds" % float(time.perf_counter() - self.start))
        self.splitter_2.addWidget(self.bioTree)
        #self.mainLogger.debug("Up to projectTree took %f seconds" % float(time.perf_counter() - self.start))
        self.splitter_2.addWidget(self.projectTree)
        self.mainLogger.debug("Adding all GUI objects took %f seconds" % float(time.perf_counter() - self.start))

        self.changeTheme()

        # Tool bar setup
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.toolBar.addWidget(spacer)
        self.toolBar.addAction(self.actionOptions)
        # Status bar setup
        self.updateUsage()
        self.statusBar().addPermanentWidget(self.memLabel)
        #self.mainLogger.debug("After StatusbarUpdate")

        # Load
        self.DEBUG()  # TODO: DELETE THIS NEPHEW

        # self.pdbWindow = displays.NGLviewer(self)

        # self.pdbWindow.show()
        del spacer

    def connectSlots(self):
        # Toolbar and MenuBar
        # FILE
        self.actionNew.triggered.connect(self.newWorkspace)
        self.actionOpen.triggered.connect(self.openWorkspace)
        self.actionImportSeq.triggered.connect(self.importSequence)
        self.actionImportAlign.triggered.connect(self.importAlignment)
        self.actionExportSeq.triggered.connect(self.exportSequence)
        self.actionExportAlign.triggered.connect(self.exportAlignment)
        self.actionSave.triggered.connect(self.saveWorkspace)
        self.actionQuit.triggered.connect(self.quit)

        # EDIT
        self.actionCopy.triggered.connect(self.copyOut)
        self.actionPaste.triggered.connect(self.pasteInto)
        # self.actionPreferences.triggered.connect(self.openPrefWindow)

        # TOOLS
        self.actionAlign.triggered.connect(self.seqDbClick)
        self.actionNewFolder.triggered.connect(self.addFolder)
        self.actionDelete.triggered.connect(self.deleteNode)

        # WINDOW
        self.actionTile.triggered.connect(self.tileWindows)
        self.actionCascade.triggered.connect(self.cascadeWindows)
        self.actionToggle_Tabs.triggered.connect(self.mdiArea.toggleTabs)  # TODO: Set as preference
        self.actionClose.triggered.connect(self.closeTab)
        self.actionClose_all.triggered.connect(self.closeAllTabs)

        # HELP
        self.actionOnThemes.triggered.connect(self.openThemeHelp)
        self.actionAbout.triggered.connect(self.showAbout)

        # Data awareness connections
        self.bioTree.doubleClicked.connect(self.seqDbClick)
        self.bioModel.dupeName.connect(self.dupeNameMsg)
        self.projectTree.doubleClicked.connect(self.alignmentDbClick)
        # self.sendParams.connect(self.optionsPane.setParams)  # This is to keep the pane in check with an opening window
        # self.optionsPane.updateParam.connect(self.nodeUpdate)  # This is to make sure the node data is up to date

        # Utility slots
        self.bioTree.generalClick.connect(self.deselectProject)
        self.bioTree.clicked.connect(self.lastClickedSeq)
        self.projectTree.generalClick.connect(self.deselectSeqs)
        self.bioModel.nameWasChanged.connect(self.postNameChange)
        self.processTimer.timeout.connect(self.updateUsage)

        # Toolbar slots
        self.actionSave_Image.triggered.connect(self.saveImage)
        self.actionOptions.toggled.connect(self.toggleOptionsPane)
        self.optionsPane.checkRuler.toggled.connect(self.toggleRuler)
        self.optionsPane.checkColors.toggled.connect(self.toggleColors)
        self.optionsPane.comboTheme.currentIndexChanged.connect(self.changeTheme)
        self.optionsPane.comboFont.currentFontChanged.connect(self.changeFont)
        self.optionsPane.spinFontSize.valueChanged.connect(self.changeFontSize)
        self.optionsPane.checkStructure.toggled.connect(self.toggleStructure)
        self.optionsPane.checkConsv.toggled.connect(self.toggleConsv)
        self.optionsPane.comboReference.currentIndexChanged.connect(self.selectReference)

        self.mdiArea.subWindowActivated.connect(self.setCurrentWindow)
        self.mdiArea.subWindowActivated.connect(self.refreshParams)

        self.optionsPane.buttonStructure.clicked.connect(self.get_UniprotId)
        self.optionsPane.checkColorDesc.toggled.connect(self.showColorDesc)

    def disconnectSlots(self):
        slots = [self.actionNew, self.actionOpen, self.actionImportSeq, self.actionImportAlign, self.actionExportSeq,
                 self.actionExportAlign, self.actionSave, self.actionQuit, self.actionCopy, self.actionPaste,
                 self.actionAlign, self.actionNewFolder, self.actionDelete, self.actionTile, self.actionCascade,
                 self.actionClose, self.actionClose_all, self.actionToggle_Tabs, self.actionAbout,
                 self.bioTree, self.bioModel, self.projectTree, self.processTimer, self.actionSave_Image,
                 self.actionOptions, self.mdiArea, #self.optionsPane.checkRuler, self.optionsPane.checkColors,
                 #self.optionsPane.comboTheme, self.optionsPane.comboFont, self.optionsPane.spinFontSize,
                 #self.optionsPane.checkConsv, self.optionsPane.comboReference,
                 self.optionsPane.buttonStructure, #self.optionsPane.checkStructure
                 ]
        for signal in slots:
            signal.disconnect()
        del slots

        #self.actionMemory.triggered.connect(self.memoryPrint)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # UTILITY METHODS
    # Represents the fundamental methods of UI interaction logic.
    def seqInit(self, seqr, folder=None):
        """
        Ran upon adding a sequence. Input is SeqRecord. Assigns a unique Window ID and adds to list of all Seqs.
        The Window ID is the core identifier for all sequence and alignment objects, and remains static for as
        long as the sequence or alignment exists.
        Creates a node in the BioModel/BioTree with the title linked to the Name
        Updating this node name only modifies the sName; the sequence ID from the FASTA remains unchanged.
        """
        wid = str(int(self.windex) + 1)
        sname = seqr.name
        # Check if title already exists, and if so, changes it.
        sname, self.titles = utilities.checkName(sname, self.titles)
        if sname != seqr.name:
            seqr.name = sname
        # Adds to the list of sequences, by its Window ID
        self.sequences[wid] = [seqr]
        node = QStandardItem(sname)
        node.setData([seqr], self.SequenceRole)
        node.setData(node.data(role=self.SequenceRole)[0].name)
        node.setData(wid, self.WindowRole)
        node.setFlags(node.flags() ^ Qt.ItemIsDropEnabled)
        if not folder:
            self.bioModel.appendRow(node)
        else:
            found = False
            for n in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
                if n.text() == folder:
                    n.appendRow(node)
                    found = True
            if not found:
                nfolder = QStandardItem(str(folder))
                nfolder.appendRow(node)
                self.bioModel.appendRow(nfolder)
                self.bioTree.setExpanded(nfolder.index(), True)
                del nfolder, n, found
        self.windex = int(wid)
        del node, wid, seqr, folder

    def seqDbClick(self):
        """
        This assesses what has been selected, and if >2 selected, adds the combo to the list of all sequences.
        Depending on whether it is a single sequence or multiple, either creates, stores and (re)displays the window or
        first makes an alignment then does all that. Ignores any folders that were included in the selection.
        Will not duplicate alignments. Creates a new window only if alignment is new.
        """
        # Items is an input dictionary for sending to clustalo
        # combo is an array of SeqRecords, sorted, to prevent creating duplicate alignments.
        self.localtime = time.perf_counter()

        items = {}
        combo = []
        # Collect the selected sequence(s)
        if self.lastClickedTree == self.bioTree:
            indices, seqs = utilities.nodeSelector(self.bioTree, self.bioModel)
            for seqr in seqs:
                items[seqr.name] = str(seqr.seq)
                combo.append(seqr)
        else:
            self.mainStatus.showMessage("Please select sequences", msecs=3000)
        combo.sort()
        aligned = self.callAlign(items)
        if items:
            if combo in self.sequences.values():
                # If an alignment with this combo has already been made...
                for key, value in self.sequences.items():
                    if combo == value:
                        # Get the window ID for this combo
                        wid = key
                        for x in range(len(combo)):
                            # Check names and rebuild window if different
                            # TODO: THIS CAN BE DELETED I THINK
                            if combo[x].name != value[x].name:
                                self.windows.pop(key)
                        try:
                            # Reopen the window, if it exists.
                            sub = self.windows[wid]
                            self.openWindow(sub)
                            del sub
                        except KeyError:
                            # Or generate a new window, if it does not.
                            sub = self.makeNewWindow(wid, aligned)
                            self.openWindow(sub)
                            del sub
            else:
                # If it hasn't been made yet, add the combo to the main list and make/open window.
                wid = str(int(self.windex) + 1)
                self.sequences[wid] = combo
                sub = self.makeNewWindow(wid, aligned)
                self.openWindow(sub)
                self.windex = self.windex + 1
                del aligned, wid, sub
        del combo, items

    def callAlign(self, seqarray):
        """
        Generates a sequence alignment using Clustal Omega in a separate thread. Currently returns with the same
        order as the nodes were clicked; need to figure out how to return with alignment order.
        This is also called upon double clicking a single sequence, but just passes it through if so.
        """
        # TODO: Do pairwise here if only 2!
        if len(list(seqarray.values())) > 1:
            # Sort the sequences to prevent duplicates and generate the alignment in a new thread.
            worker = utilities.AlignThread(self, seqarray, seqtype=3, num_threads=self.threadpool.maxThreadCount())
            worker.start()
            worker.finished.connect(worker.deleteLater)
            worker.finished.connect(worker.quit)
            worker.wait()
            aligned = worker.aligned
            del worker
        #  elif len(list(seqarray.values())) == 2:
        else:
            aligned = seqarray  # send single sequence
        return aligned

    def alignmentDbClick(self):
        """ Simple method that confirms you didn't click on a folder, then opens the window """
        # Only a single item is selectable at once in the alignment tree.
        item = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
        if item.data(role=self.WindowRole):
            sub = self.windows[item.data(role=self.WindowRole)]
            sub.setWindowTitle(item.data())
            self.openWindow(sub)
            del sub
        del item

    def makeNewWindow(self, wid, ali, nonode=False):
        """
        Method for generating a new subwindow in the MDI area. Usually called upon double-clicking
        a sequence node, or multiple sequence nodes, if those are not already in the form of a window.
        After generation, stores the window in the list of all windows by its window ID.
        """
        sub = widgets.MDISubWindow(wid)
        widget = AlignSubWindow(ali, self.optionsPane.params)
        sub.setWidget(widget)
        if len(ali.keys()) > 1:
            # If alignment:
            if nonode:
                self.windows[wid] = sub
            else:
                seqrs = []
                for key, value in ali.items():
                    for stored_seq in self.sequences.values():
                        if str(value).replace('-', "") == str(stored_seq[0].seq):
                            # print("USING ORIGINAL ID", stored_seq[0].id)
                            sid = stored_seq[0].id
                            break
                        else:
                            # print("USING NAME", key)
                            sid = key
                    # print("Adding sequence to alignment: ", key, sid)
                    del stored_seq
                    seqr = models.SeqR(Seq(value, generic_protein), name=key, id=sid)
                    seqrs.append(seqr)
                aliR = MultipleSeqAlignment(seqrs)
                node = QStandardItem(sub.windowTitle())
                node.setData(sub.windowTitle())
                node.setData(aliR, self.SequenceRole)
                node.setData(wid, self.WindowRole)
                # node.setData(self.params.copy(), role=self.AlignDisplayRole)
                self.projectRoot.appendRow(node)
                self.projectTree.setExpanded(node.parent().index(), True)
                self.windows[wid] = sub
                self.projectModel.updateWindows(self.windows)
                del node
                del aliR
                del seqr
                del seqrs
        else:
            # If sequence:
            seq = self.sequences[wid][0]
            sub.setWindowTitle(seq.name)
            self.windows[wid] = sub
            self.bioModel.updateNames(self.titles)
            self.bioModel.updateWindows(self.windows)
            del seq
            del widget
        return sub

    def openWindow(self, sub):
        """
        Checks to see if a window is open already.
        If it is not, reopens the window. If it is, gives focus.
        """
        self.localtime = time.perf_counter()
        found = False
        for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            if node.data(self.WindowRole) == sub.wid:
                if node.data(self.StructureRole):
                    sub.widget().addStructure(node.data(self.StructureRole), node.data(self.SequenceRole)[0])
                    found = True
        if not found:
            seq = None
            for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
                if node.data(self.WindowRole) == sub.wid:
                    if node.data(self.StructureRole):
                        seq = node.data(self.SequenceRole)[0]
                        sub.widget().addStructure(node.data(self.StructureRole), seq)
                        found = True
                        del seq
        del found, node
        sub.widget().setParams(self.optionsPane.params)
        if sub.mdiArea() != self.mdiArea:
            self.mainLogger.debug("Adding window to MDI Area; creation took %f seconds" %
                                  float(time.perf_counter() - self.localtime))
            # print(sub.widget().params)
            # self.sendParams.emit(sub.widget().params.copy())
            self.mdiArea.addSubWindow(sub)
            #print("SUBWINDOWS: %s" % len(self.mdiArea.subWindowList()))
        else:
            sub.show()
            self.mdiArea.setActiveSubWindow(sub)
        del sub

    def maybeClose(self):
        """
        Simple confirmation upon making a new workspace or closing the program.
        Currently not using as the behavior is inconsistent with the window 'X' button.
        """
        qDialog = QuitDialog(self)
        confirm = qDialog.exec()
        if confirm == 1:
            result = self.saveWorkspace()
            if result:
                return True
        elif confirm == 2:
            return False
        else:
            return None

    def quit(self):
        # TODO: Make this consistent and more reliable.
        # confirm = self.maybeClose()
        # if confirm is not None:
        #    self.close()
        self.close()

    def lastClickedSeq(self):
        """Records the last clicked sequence, for copying/pasting etc"""
        self.bioModel.updateNames(self.titles)
        self.bioModel.updateLastClicked(self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0]))

    def deselectProject(self):
        """Function to unselect the opposite tree when clicked. Prevents confusion when pasting"""
        self.projectTree.clearSelection()
        self.lastClickedTree = self.bioTree

    def deselectSeqs(self):
        """Function to unselect the opposite tree when clicked. Prevents confusion when pasting"""
        self.bioTree.clearSelection()
        self.lastClickedTree = self.projectTree

    def preNameChange(self):
        self.pruneNames()

    def postNameChange(self, name):
        """
        When a sequence name is changed, I have to clean up the title list,
        update that list in the bioModel, and update all the titles in the alignment window.
        """
        print("Post name change to %s" % name)
        self.pruneNames()
        self.bioModel.updateNames(self.titles)

        # TODO: This is the best way to potentially rename alignments.
        # for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
        #    wid = node.data(role=self.WindowRole)
        #    if wid:
        #        sub = self.windows[wid]
        #        subseqs = sub.seqs()
        #        seqrs = node.data(role=self.SequenceRole)
        #        for seqr in seqrs:
        #            if seqr.name

    def pruneNames(self):
        """ This checks which names are in "TITLES" and deletes if they were not. """
        names = []
        pruned = []
        for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            if node.data(role=self.SequenceRole):
                names.append(node.text())
                for key, value in self.sequences.items():
                    if node.data(role=self.SequenceRole) == value:
                        value[0].name = node.data(role=self.SequenceRole)[0].name
        for title in self.titles:
            if title not in names:
                pruned.append(title)
        self.mainLogger.debug("Tree names: " + str(names) + " vs. Stored names: " + str(self.titles))
        self.titles = [x for x in self.titles and names if x not in pruned]
        self.mainLogger.debug("Removed " + str(pruned) + ", leaving " + str(self.titles))
        del names, pruned, title, node

    def updateUsage(self):
        """ Simple method that updates the status bar process usage statistics on timer countdown """
        # self.tr.print_diff()
        mem = self.mainProcess.memory_full_info().uss / 1000000
        cpu = self.mainProcess.cpu_percent()
        self.memLabel.setText("CPU: " + str(cpu) + " % | RAM: " + str(round(mem, 2)) + " MB")


    #def memoryPrint(self):
     #   print("MEMORY OUTPUT")
      #  #cb = refbrowser.ConsoleBrowser(self._currentWindow.widget(), maxdepth=5, str_func=utilities.output_function)
       # self.tr = tracker.SummaryTracker()
        #self.tr.print_diff()
        #cb.print_tree()
        #sum2 = summary.summarize(muppy.get_objects())
        #diff = summary.get_diff(self.sum, sum2)
        #summary.print_(diff)
        #self.sum = sum2

        #objgraph.show_growth()
        #objgraph.show_chain(
        #    objgraph.find_backref_chain(
        #        random.choice(objgraph.by_type(dict)),
        #        objgraph.is_proper_module), )

    def newWorkspace(self):
        """ Open a new workspace; ask first if you want to save what you are working on. """
        result = self.maybeClose()
        if result is not None:
            self.mainLogger.info("CLEARING WORKSPACE")
            #self.deleteLater()
            #newWindow = Linnaeo()
            #newWindow.setWindowTitle("Linnaeo")
            #icon = QIcon(":/icons/linnaeo_full.ico")
            #newWindow.setWindowIcon(icon)
            #qApp._window = newWindow
            #newWindow.show()
            self.mdiArea.closeAllSubWindows()
            self.start = time.perf_counter()
            self.disconnectSlots()
            self.guiSet()
            self.connectSlots()
            #self.destroy()
            del result
        else:
            del result

    def openWorkspace(self):
        """ Function for opening a saved workspace. Goes through and reloads in reverse of how it is saved. """
        sel = QFileDialog.getOpenFileName(parent=self, caption="Open Workspace", directory=QDir.homePath(),
                                          filter="Linnaeo Workspace (*.lno);;Any (*)")
        if sel != ('', ''):
            filename = sel[0]
            self.mainLogger.debug("Opening file: " + str(filename))
            file = QFile(filename)
            file.open(QIODevice.ReadOnly)
            fstream = QDataStream(file)
            self.mainLogger.info("Starting restore")
            self.mainStatus.showMessage("Rebuilding workspace, please wait...")
            sequences = fstream.readQVariantHash()
            titles = fstream.readQVariantList()
            windex = fstream.readUInt32()
            windows = {}
            self.mdiArea.closeAllSubWindows()
            newBModel = widgets.ItemModel(windows, seqTree=True)
            newPModel = widgets.ItemModel(windows)
            self.restore_tree(newBModel.invisibleRootItem(), fstream)
            self.restore_tree(newPModel.invisibleRootItem(), fstream)
            self.mainStatus.showMessage("Rebuild complete!")
            self.mdiArea.closeAllSubWindows()
            self.start = time.perf_counter()
            self.disconnectSlots()
            self.guiSet(trees=[newBModel, newPModel], data=[sequences, titles, windex])
            self.connectSlots()
            del sel, filename, file, fstream, sequences, titles, windex, windows, newBModel, newPModel
            #self.destroy()
        else:
            del sel


class LinnaeoApp(QApplication):
    """ Custom QApplication that I made to detect events. Currently not needed."""

    def __init__(self, *args):
        super().__init__(*args)
        """
        # SAVED THIS FOR LATER!
        # Read in config (linux)
        config = configparser.ConfigParser()
        try:
            config.read(str(Path.home())+"/.linnaeo/config.ini")
            # Open last used workspace automatically.
            if config['RECENTS']['LAST'] != "":
                last = config['RECENTS']['LAST']
                f = QFile(last)
                f.open(QIODevice.ReadOnly)
                model = workspace.WorkspaceModel()
                f.close()
                self.workspaceTree.setModel(model)
        except:
            print("No config file found!")
            """
        self.mode = "Light"
        self.fonts = QFontDatabase()
        self.defFontId = self.fonts.addApplicationFont(':/fonts/Default-Noto.ttf')
        self.defFontId2 = self.fonts.addApplicationFont(':/fonts/LiberationMono.ttf')

        # For some reason, Noto does not work on Mac or Windows. But does for SS display???
        self.defFont = QFont(self.fonts.applicationFontFamilies(self.defFontId2)[0], 10) \
            if sys.platform in ['win32', 'darwin'] else QFont(self.fonts.applicationFontFamilies(self.defFontId)[0], 10)

    def setMode(self, mode):
        self.mode = mode
