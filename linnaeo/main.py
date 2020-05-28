#!/usr/bin/env python3

# Bioscience components
import logging
import os
import time

import psutil
# PyQt components
from PyQt5.QtCore import Qt, QThreadPool, pyqtSignal, QCoreApplication
from PyQt5.QtGui import QStandardItem
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QAbstractItemView

# Internal components
import linnaeo
from linnaeo.resources import linnaeo_rc
from linnaeo.classes import widgets, utilities, methods
from linnaeo.classes.displays import QuitDialog, AlignSubWindow
from linnaeo.ui import linnaeo_ui


class Linnaeo(QMainWindow, methods.Slots, methods.Debug, linnaeo_ui.Ui_MainWindow):
    edgeClick = pyqtSignal()
    """
    Constructor for the Main Window of the Sherlock App
    Contains all the user interface functions, as well as underlying code for major features.
    The bulk of the program is located here.
    Please see classes.methods for additional methods for this class. 
    """

    def __init__(self, *args, **kwargs):

        super(self.__class__, self).__init__(*args, **kwargs)
        self.start = linnaeo.start_time
        # Initialize UI
        self.beingClicked = True
        self.setAttribute(Qt.WA_QuitOnClose)
        self.setupUi(self)  # Built by PyUic5 from my main window UI file

        # System instants; status bar and threads
        self.memLabel = QLabel()
        self.mainProcess = psutil.Process(os.getpid())
        # self.processTimer = QTimer()
        self.processTimer = utilities.TimerThread()
        self.mainLogger = logging.getLogger("Main")
        self.threadpool = QThreadPool()
        self.mainLogger.info("Threading with a maximum of %d threads" % self.threadpool.maxThreadCount())

        # Project instants and inherent variables for logic.
        self.localtime = 0
        self.SequenceRole = Qt.UserRole + 2
        self.WindowRole = Qt.UserRole + 3
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

        # MDI Window
        self.mdiArea = widgets.MDIArea()
        self.gridLayout_2.addWidget(self.mdiArea)

        # Tree stuff
        self.bioTree = widgets.TreeView()
        self.projectTree = widgets.TreeView()
        self.mainLogger.debug("Finished initializing, took %f seconds" % float(time.perf_counter() - self.start))

        # Other functions
        self.guiSet()
        self.mainLogger.debug("Gui initializing complete, took %f seconds" % float(time.perf_counter() - self.start))
        self.guiFinalize()
        self.mainLogger.debug("Finalized gui, took %f seconds" % float(time.perf_counter() - self.start))
        self.connectSlots()
        self.mainLogger.debug("Slots connected after %f seconds" % float(time.perf_counter() - self.start))
        self.mainLogger.debug("Setup took %f seconds" % float(time.perf_counter() - self.start))
        self.setMouseTracking(True)

    def guiSet(self, trees=None, data=None):
        """ Initialize GUI with default parameters. """
        self.lastClickedTree = None
        self.lastAlignment = {}
        self.windows = {}  # Windows stored as { windex : MDISubWindow }
        self.windex = 0  # Acts as identifier for tracking alignments (max 2.1 billion)
        self.sequences = {}  # Stored as {WindowID:[SeqRecord(s)] }
        self.titles = []  # maintains a list of sequence titles to confirm uniqueness
        self.mainLogger.debug("guiSet took took %f seconds" % float(time.perf_counter() - self.start))

        if trees:
            # For loading a window on File>Open
            self.mainLogger.info("Loading saved workspace!")
            self.bioModel, self.projectModel = trees
            self.sequences, self.titles, self.windex = data
            self.windex = int(self.windex)
            self.bioRoot = self.bioModel.invisibleRootItem().child(0)  # TODO: ELIMINATE USE OF ROOT. Use InvsRoot
            self.projectRoot = self.projectModel.invisibleRootItem().child(0)

            self.windows = {}
            self.rebuildTrees()
            # TODO: Figure out why tree is not expanding!
            #print("SEQS: ", self.sequences)
            #print("TITLES: ", self.titles)
            #print("WINDEX: ", self.windex)
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
        self.installEventFilter(self)
        self.statusBar().setSizeGripEnabled(False)

        #self.splitter_2.splitterMoved.connect(self.setSizing)

    def eventFilter(self, obj, event):
        """ Designed to capture the edge mouse click upon resizing """
        #print(event, event.type())
        if event.type() in [99]:
            #print("Detected edge")
            self.edgeClick.emit()
        return super().eventFilter(obj, event)

    def resizeEvent(self, event):
    #    # TODO: THIS IS THE ONLY ONE THAT DOESN'T TURN OFF COLOR WHEN RESIZING
        # Just flashes of and on...
        # self.setSizing()
        return super().resizeEvent(event)

    def guiFinalize(self):
        # Tree setup
        self.mainLogger.debug("Up to selectionMode took %f seconds" % float(time.perf_counter() - self.start))
        self.bioTree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.mainLogger.debug("Up to adding bioTree took %f seconds" % float(time.perf_counter() - self.start))
        self.splitter_2.addWidget(self.bioTree)
        self.mainLogger.debug("Up to projectTree took %f seconds" % float(time.perf_counter() - self.start))
        self.splitter_2.addWidget(self.projectTree)
        self.mainLogger.debug("Adding all GUI objects took %f seconds" % float(time.perf_counter() - self.start))

        # Status bar setup
        self.updateUsage()
        self.statusBar().addPermanentWidget(self.memLabel)
        self.mainLogger.debug("After StatusbarUpdate")
        # self.processTimer.setInterval(1000)
        # self.processTimer.start()

        self.DEBUG()

    def connectSlots(self):
        # Toolbar and MenuBar
        # FILE
        self.actionNew.triggered.connect(self.newWorkspace)
        self.actionOpen.triggered.connect(self.openWorkspace)
        self.actionAdd.triggered.connect(self.newSequence)
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
        self.actionAbout.triggered.connect(self.showAbout)

        # Data awareness connections
        self.bioTree.doubleClicked.connect(self.seqDbClick)
        self.bioModel.dupeName.connect(self.dupeNameMsg)
        self.projectTree.doubleClicked.connect(self.alignmentDbClick)

        # Utility slots
        self.bioTree.generalClick.connect(self.deselectProject)
        self.bioTree.clicked.connect(self.lastClickedSeq)
        self.projectTree.generalClick.connect(self.deselectSeqs)
        # self.bioModel.nameChanging.connect(self.preNameChange)
        self.bioModel.nameChanged.connect(self.postNameChange)
        self.processTimer.timeout.connect(self.updateUsage)

        self.actionRulers.triggered.connect(self.toggleRulers)
        self.actionColors.triggered.connect(self.toggleColors)
        self.actionSave_Image.triggered.connect(self.saveImage)
        LinnaeoApp.instance().barClick.connect(self.setSizing)
        self.edgeClick.connect(self.setSizing)
        self.mdiArea.subWindowActivated.connect(self.setCurrentWindow)
        self.actionBigger.triggered.connect(self.increaseTextSize)
        self.actionSmaller.triggered.connect(self.decreaseTextSize)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # UTILITY METHODS
    # Represents the fundamental methods of UI interaction logic.
    def seqInit(self, seqr):
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
        self.bioModel.appendRow(node)
        self.windex = int(wid)

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
        self.mainLogger.debug("Beginning window creation at %s" % float(self.localtime))
        items = {}
        combo = []
        # Collect the selected sequence(s)
        for index in self.bioTree.selectedIndexes():
            # Need to make a dictionary of { Name:Sequence } for sending to ClustalO.
            # Only does the thing if there is a sequence present in the node.
            if self.bioModel.itemFromIndex(index).data(role=self.SequenceRole):
                seqr = self.bioModel.itemFromIndex(index).data(role=self.SequenceRole)[0]
                items[seqr.name] = str(seqr.seq)
                combo.append(seqr)
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
                        except KeyError:
                            # Or generate a new window, if it does not.
                            sub = self.makeNewWindow(wid, aligned)
                            self.openWindow(sub)
            else:
                # If it hasn't been made yet, add the combo to the main list and make/open window.
                wid = str(int(self.windex) + 1)
                self.sequences[wid] = combo
                sub = self.makeNewWindow(wid, aligned)
                self.openWindow(sub)
                self.windex = self.windex + 1

    def callAlign(self, seqarray):
        """
        Generates a sequence alignment using Clustal Omega in a separate thread. Currently returns with the same
        order as the nodes were clicked; need to figure out how to return with alignment order.
        This is also called upon double clicking a single sequence, but just passes it through if so.
        """
        # TODO: Do pairwise here if only 2!
        if len(list(seqarray.values())) > 1:
            # Sort the sequences to prevent duplicates and generate the alignment in a new thread.
            worker = utilities.AlignThread(seqarray, seqtype=3, num_threads=self.threadpool.maxThreadCount())
            worker.start()
            worker.wait()
            aligned = worker.aligned
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
            self.openWindow(sub)

    def makeNewWindow(self, wid, ali, nonode=False):
        """
        Method for generating a new subwindow in the MDI area. Usually called upon double-clicking
        a sequence node, or multiple sequence nodes, if those are not already in the form of a window.
        After generation, stores the window in the list of all windows by its window ID.
        """
        sub = widgets.MDISubWindow()
        widget = AlignSubWindow(ali)
        sub.setWidget(widget)
        if len(ali.keys()) > 1:
            if nonode:
                self.windows[wid] = sub
            else:
                node = QStandardItem(sub.windowTitle())
                node.setData(sub.windowTitle())
                node.setData(self.sequences[wid], self.SequenceRole)
                node.setData(wid, self.WindowRole)
                self.projectRoot.appendRow(node)
                self.projectTree.setExpanded(node.parent().index(), True)
                self.windows[wid] = sub
                self.projectModel.updateWindows(self.windows)
        else:
            seq = self.sequences[wid][0]
            sub.setWindowTitle(seq.name)
            self.windows[wid] = sub
            self.bioModel.updateNames(self.titles)
            self.bioModel.updateWindows(self.windows)
        return sub

    def openWindow(self, sub):
        """
        Checks to see if a window is open already.
        If it is not, reopens the window. If it is, gives focus.
        """
        if sub.mdiArea() != self.mdiArea:
            self.mainLogger.debug("Adding window to MDI Area; \
                                creation took %f seconds" % float(time.perf_counter() - self.localtime))
            self.mdiArea.addSubWindow(sub)
        else:
            sub.show()
            self.mdiArea.setActiveSubWindow(sub)

    def rebuildTrees(self):
        """
        Run upon loading a saved file.
        At this point, the sequences and the alignments both have Window IDs applied -- but the windows
        no longer exist. This rebuilds the windows using the saved window IDs and updates the respective models.
        """
        for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            if node.data(role=self.SequenceRole):
                self.mainLogger.info("Loading sequence: "+node.data())
                ali = {}  # empty dict needed to send to open window
                wid = node.data(role=self.WindowRole)
                seqr = node.data(role=self.SequenceRole)[0]
                ali[seqr.name] = str(seqr.seq)
                self.makeNewWindow(wid, ali, nonode=True)
                #self.bioTree.setExpanded(self.bioModel.indexFromItem(node), True)

        for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
            if node.data(role=self.SequenceRole):
                self.mainLogger.info("Loading alignment: "+node.data())
                seqs = {}
                wid = node.data(role=self.WindowRole)
                seqr = node.data(role=self.SequenceRole)
                for seq in seqr:
                    seqs[seq.name] = str(seq.seq)
                worker = utilities.AlignThread(seqs, seqtype=3, num_threads=self.threadpool.maxThreadCount())
                worker.start()
                worker.wait()
                ali = worker.aligned
                self.makeNewWindow(wid, ali, nonode=True)
                #self.projectTree.setExpanded(self.projectModel.indexFromItem(node), True)
        self.bioModel.updateWindows(self.windows)
        self.projectModel.updateWindows(self.windows)
        self.mainLogger.debug("Regenerating windows took took %f seconds" % float(time.perf_counter() - self.start))

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
        #confirm = self.maybeClose()
        #if confirm is not None:
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

    def postNameChange(self):
        """
        When a sequence name is changed, I have to clean up the title list,
        update that list in the bioModel, and update all the titles in the alignment window.
        """
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
        """
        CONNECTS TO SIGNAL
        This checks which names are in "TITLES" and deletes if they were not.
        """
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

    def updateUsage(self):
        """ Simple method that updates the status bar process usage statistics on timer countdown """
        mem = self.mainProcess.memory_info().rss / 1000000
        cpu = self.mainProcess.cpu_percent()
        self.memLabel.setText("CPU: " + str(cpu) + " % | RAM: " + str(round(mem, 2)) + " MB")


class LinnaeoApp(QApplication):
    """ Custom QApplication that detects when the border is being pressed. """
    barClick = pyqtSignal()

    def __init__(self, *args):
        super().__init__(*args)
        self.focusChanged.connect(self.focusChange)

        QCoreApplication.setAttribute(Qt.AA_CompressHighFrequencyEvents)
        self.installEventFilter(self)
        self.last = None
        #self.barClick.connect(self.setSizing)

    def focusChange(self):
        print("FOCUS CHANGED")

    def eventFilter(self, obj, event):
        if event.type() in [174,175]:
            if event.type() != self.last:
                print("Unique event!", event.type())
                print("Last observed to be: ", self.last)
                self.barClick.emit()
                self.last = event.type()
        return super().eventFilter(obj, event)

    """
     def event(self, event):
         if event.type() in [2,3]:
             print("EVENT FROM QAPP")
         #print(event, event.type())
         return super().event(event)
     """
