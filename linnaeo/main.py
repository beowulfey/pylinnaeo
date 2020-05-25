#!/usr/bin/env python3

# Bioscience components
import logging
import os
import time

import Bio.Seq as Bseq
import psutil
from Bio import SeqIO
from Bio.Alphabet import generic_protein

# PyQt components
from PyQt5.QtCore import Qt, QThreadPool, QDir, QFile, QIODevice, QDataStream, pyqtSignal
from PyQt5.QtGui import QStandardItem, QCursor
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QAbstractItemView, QFileDialog

# Internal components
import linnaeo
from linnaeo.classes import models, widgets, utilities
from linnaeo.classes.displays import QuitDialog, AlignSubWindow, AboutDialog
from linnaeo.ui import linnaeo_ui


class Linnaeo(QMainWindow, linnaeo_ui.Ui_MainWindow):
    edgeClick = pyqtSignal()
    """
    Constructor for the Main Window of the Sherlock App
    Contains all the user interface functions, as well as underlying code for major features.
    The bulk of the program is located here.
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
            print("SEQS: ", self.sequences)
            print("TITLES: ", self.titles)
            print("WINDEX: ", self.windex)
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
        self.bioTree.setExpanded(self.bioModel.invisibleRootItem().index(), True)
        self.bioModel.setHorizontalHeaderLabels(["Sequences"])
        self.projectModel.setHorizontalHeaderLabels(["Alignments"])
        self.projectTree.setExpanded(self.projectModel.invisibleRootItem().index(), True)
        self.installEventFilter(self)

    def event(self, event):
        #print(event, event.type())
        return super().event(event)

    def eventFilter(self, obj, event):
        """ Designed to capture the edge mouse click upon resizing """
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

    def increaseTextSize(self):
        self._currentWindow.widget().increaseFont()
        self.mainStatus.showMessage("Setting font size to %s " % self._currentWindow.widget().getFontSize(), msecs=1000)

    def decreaseTextSize(self):
        self._currentWindow.widget().decreaseFont()
        self.mainStatus.showMessage("Setting font size to %s " % self._currentWindow.widget().getFontSize(), msecs=1000)

    def setCurrentWindow(self):
        self._currentWindow = self.mdiArea.activeSubWindow()

    def showAbout(self):
        qDialog = AboutDialog(self)
        qDialog.exec()
        """if confirm == 1:
            result = self.saveWorkspace()
            if result:
                print("YES")
                return True
        elif confirm == 2:
            print("NO")
            return False
        else:
            return None
            """

    def setSizing(self):
        if not self.beingClicked:
            self.beingClicked = True
            if self._currentWindow and self._currentWindow.isMaximized():  # and self.mdiArea.activeSubWindow().isMaximized():
                # print("REDRAWING FRAME FROM MAIN")
                self._currentWindow.widget().userIsResizing = True
                self._currentWindow.widget().seqArrange(color=False, rulers=False)
        elif self.beingClicked:
            self.beingClicked = False
            if self._currentWindow and self._currentWindow.isMaximized():
                # print("DONE REDRAWING FROM MAIN")
                self._currentWindow.widget().userIsResizing = False
                self._currentWindow.widget().seqArrange()

    def saveImage(self):
        file = QFileDialog.getSaveFileName(self, "Save as...", "name",
                                            "PNG (*.png);; BMP (*.bmp);;TIFF (*.tiff *.tif);; JPEG (*.jpg *.jpeg)");
        self._currentWindow.widget().seqArrange(color=True, rulers=True)
        self._currentWindow.widget().grab().save(file[0]+file[1][-5:-1])

    def toggleRulers(self):
        self._currentWindow.widget().toggleRulers()

    def toggleColors(self):
        self._currentWindow.widget().toggleColors()

    # SINGLE-USE SLOTS
    def newWorkspace(self):
        result = self.maybeClose()
        if result is not None:
            self.mainLogger.info("CLEARING WORKSPACE")
            self.mdiArea.closeAllSubWindows()
            self.guiSet()

    def restore_tree(self, parent, datastream, num_childs=None):
        if not num_childs:
            #print("First line: reading UInt32")
            num_childs = datastream.readUInt32()
            #print(num_childs)
        for i in range(0, num_childs):
            #print("Reading node")
            child = QStandardItem()
            child.read(datastream)
            parent.appendRow(child)
            #print(child.data(role=Qt.UserRole + 1))
            #print(child.data(role=Qt.UserRole + 2))
            #print(child.data(role=Qt.UserRole + 3))
            num_childs = datastream.readUInt32()
            if num_childs > 0:
                #print("reading children")
                self.restore_tree(child, datastream, num_childs)

    def openWorkspace(self):
        sel = QFileDialog.getOpenFileName(parent=self, caption="Open Workspace", directory=QDir.homePath(),
                                          filter="Linnaeo Workspace (*.lno);;Any (*)")
        filename = sel[0]
        self.mainLogger.debug("Opening file: " + str(filename))
        file = QFile(filename)
        file.open(QIODevice.ReadOnly)
        fstream = QDataStream(file)
        self.mainLogger.info("Starting restore")
        sequences = fstream.readQVariantHash()
        print("Sequences: ", sequences)
        titles = fstream.readQVariantList()
        print("Titles: ", titles)
        windex = fstream.readUInt32()
        print("Windex: ", windex)
        windows = {}
        self.mdiArea.closeAllSubWindows()
        newBModel = widgets.ItemModel(windows, seqTree=True)
        newPModel = widgets.ItemModel(windows)
        self.restore_tree(newBModel.invisibleRootItem(), fstream)
        self.restore_tree(newPModel.invisibleRootItem(), fstream)
        self.guiSet(trees=[newBModel, newPModel], data=[sequences, titles, windex])

    def newSequence(self):
        sel = QFileDialog.getOpenFileName(parent=self, caption="Load a Single Sequence", directory=QDir.homePath(),
                                          filter="Fasta File (*.fa *.fasta);;Any (*)")
        seqr = None
        filename = sel[0]
        badfile = True
        if filename[-3:] == ".fa":
            badfile = False
        elif filename[-6:] == ".fasta":
            badfile = False
        if not badfile:
            try:
                seq = SeqIO.read(filename, "fasta")
                if seq:
                    #print("name: ", seq.name)
                    #print("id: ", seq.id)
                    #print("desc: ", seq.description)
                    seqr = models.SeqR(seq.seq, name=seq.name, id=seq.id, description=seq.description)
                    self.seqInit(seqr)
            except:
                self.mainStatus.showMessage("ERROR -- Please check file", msecs=3000)
        else:
            self.mainStatus.showMessage("Please only add fasta file!", msecs=1000)

    def saveWorkspace(self):
        sel = QFileDialog.getSaveFileName(parent=self, caption="Save Workspace", directory=QDir.homePath(),
                                          filter="Linnaeo Workspace (*.lno);;Any (*)")
        filename = sel[0]
        if filename:
            if filename[-4:] != ".lno":
                filename = str(filename + ".lno")
            file = QFile(filename)
            file.open(QIODevice.WriteOnly)
            out = QDataStream(file)
            self.mainLogger.debug("Beginning file save")
            #print("TREE DATA")
            #self.queryTrees()
            #print("Writing Sequences, Titles, Windex")
            #print(self.sequences)
            #print(self.titles)
            #print(type(self.windex), self.windex)
            out.writeQVariantHash(self.sequences)
            out.writeQVariantList(self.titles)
            out.writeUInt32(self.windex)
            #####################################
            # Sequence View
            ###
            print("SEQUENCE TREE")
            print("Invs.Root Children: ", self.bioModel.invisibleRootItem().rowCount())
            out.writeUInt32(self.bioModel.invisibleRootItem().rowCount())
            for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
                print("Text: ", str(node.text()))
                print("Data1: ", str(node.data()))
                print("Data2: ", str(node.data(role=self.SequenceRole)))
                print("Data3: ", str(node.data(role=self.WindowRole)))
                node.write(out)
                out.writeUInt32(node.rowCount())
            #####################################
            # Alignment View
            # Does not save any metadata! Only the two sequences
            # So I can't save window options at the moment.
            # TODO: Consider adding "window modes" to the node.
            print("ALIGNMENT TREE")
            print("Invs.Root Children: ", self.bioModel.invisibleRootItem().rowCount())
            out.writeUInt32(self.projectModel.invisibleRootItem().rowCount())
            for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
                print("Text: ", str(node.text()))
                print("Data1: ", str(node.data()))
                print("Data2: ", str(node.data(role=self.SequenceRole)))
                print("Data3: ", str(node.data(role=self.WindowRole)))
                node.write(out)
                out.writeUInt32(node.rowCount())
            self.mainLogger.debug("Save complete")
            file.flush()
            file.close()
            return True
        else:
            #print("No filename chosen; canceling")
            return False

    def copyOut(self):
        if self.lastClickedTree == self.bioTree:
            seqs = []
            indices = self.lastClickedTree.selectedIndexes()
            for index in indices:
                node = self.bioModel.itemFromIndex(index)
                seqr = node.data(role=self.SequenceRole)[0]
                print(seqr.format("fasta"))
                seqs.append(seqr.format("fasta"))
            QApplication.clipboard().setText("".join(seqs))
            if len(seqs) > 1:
                self.mainStatus.showMessage("Copied sequences to clipboard!", msecs=1000)
            elif len(seqs) == 1:
                self.mainStatus.showMessage("Copied sequence to clipboard!", msecs=1000)

        elif self.lastClickedTree == self.projectTree:
            self.mainStatus.showMessage("Please use File>Export to save alignments!", msecs=1000)

    def pasteInto(self):
        """
        Only accepts single FASTA paste right now
        """
        # FASTA DETECTION
        if self.lastClickedTree == self.bioTree:
            nline = None
            sid = None
            desc = None
            name = None
            bars = []
            spaces = []
            try:
                clip = QApplication.clipboard().text()
                if clip[0] == ">":
                    for index in range(len(clip)):
                        if not nline and clip[index] == "\n":
                            nline = clip[:index]
                            seq = clip[index + 1:]
                    seq = seq.replace('\n', '')
                    seq = seq.replace('\r', '')
                    bseq = Bseq.MutableSeq(seq)

                    for index in range(len(nline)):
                        if nline[index] == " ":
                            spaces.append(index)
                        if nline[index] == "|":
                            bars.append(index)
                    if len(bars) > 0:
                        sid = nline[bars[0] + 1:bars[1]]
                        self.mainLogger.debug("Extracted ID for sequence: ", sid)
                        desc = nline[bars[1] + 1:]
                    elif not sid:
                        sid = nline[1:9]
                        name = sid

                    if sid and bseq:
                        seqr = models.SeqR(bseq, id=sid, name=sid)
                        if desc:
                            seqr.description = desc
                        self.seqInit(seqr)
                else:
                    self.mainStatus.showMessage("Please only paste in FASTA format!", msecs=1000)
            except:
                self.mainStatus.showMessage("Please only paste in FASTA format!", msecs=1000)
        else:
            self.mainStatus.showMessage("Please choose paste destination", msecs=1000)

    def addFolder(self):
        new = QStandardItem("New Folder")
        new.setData("New Folder")
        try:
            node = self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0])
            if not node.data(role=self.WindowRole):
                node.appendRow(new)
            else:
                self.bioModel.appendRow(new)
        except IndexError:
            try:
                node = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
                if not node.data(role=self.WindowRole):
                    node.appendRow(new)
                else:
                    self.projectModel.appendRow(new)
            except IndexError:
                if self.lastClickedTree:
                    self.mainLogger.debug("Nothing selected so adding folder at last click")
                    self.lastClickedTree.model().appendRow(new)
                else:
                    self.mainStatus.showMessage("Select a location first!", msecs=1000)

    def deleteNode(self):
        for index in self.lastClickedTree.selectedIndexes():
            node = self.lastClickedTree.model().itemFromIndex(index)
            wid = node.data(role=self.WindowRole)
            try:
                sub = self.windows[wid]
                if sub:
                    self.mainLogger.debug("Deleting node from tree: " + str(sub.windowTitle()))
                    sub.close()
                    self.windows.pop(wid)
                    self.sequences.pop(wid)
                    self.pruneNames()
                    self.bioModel.updateNames(self.titles)
                    self.bioModel.updateWindows(self.windows)
                    self.projectModel.updateNames(self.titles)
                    self.projectModel.updateWindows(self.windows)
            except KeyError:
                pass
            except ValueError:
                pass
            self.lastClickedTree.model().removeRow(index.row(), index.parent())
            if self.mdiArea.tabbed:
                self.mdiArea.activeSubWindow().showMaximized()

    def tileWindows(self):
        if self.mdiArea.tabbed:
            self.mdiArea.toggleTabs()
        self.mdiArea.tileSubWindows()

    def cascadeWindows(self):
        if self.mdiArea.tabbed:
            self.mdiArea.toggleTabs()
        self.mdiArea.cascadeSubWindows()

    def closeTab(self):
        """ This duplicates inherent functionality of QMDISubWindow, so I can't add shortcut to menu"""
        try:
            self.mdiArea.activeSubWindow().close()
        except AttributeError:
            pass

    def closeAllTabs(self):
        self.mdiArea.closeAllSubWindows()

    def callAlign(self, seqarray):
        # Process the sequences, and generate an alignment if >2
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

    # DATA AWARENESS SLOTS
    def seqDbClick(self):
        """
        This assesses what has been selected, adds that to the parent list of sequences,
        and depending on whether it is a single sequence or multiple,
        either creates, stores and (re)displays the window or first makes an alignment then does all that.
        Ignores any folders that were included in the selection.
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

    def alignmentDbClick(self):
        # Checks if not a folder first, then:
        # Gets the selected item (only single selection allowed), and opens the window
        item = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
        if item.data(role=self.WindowRole):
            sub = self.windows[item.data(role=self.WindowRole)]
            self.openWindow(sub)

    def dupeNameMsg(self):
        self.mainStatus.showMessage("Please choose a unique name!", msecs=1000)

    # MAIN METHODS
    def seqInit(self, seqr):
        """
        Input is SeqRecord.
        Assigns a unique WID and adds to list of all Seqs.
        Creates a node with the title linked to the sName and adds to the Sequence Tree
        Updating this node name should modify the sName too.
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

    def makeNewWindow(self, wid, ali, nonode=False):
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
        Also refreshes the title.
        """
        if sub.mdiArea() != self.mdiArea:
            self.mainLogger.debug("Finished window at %f" % time.perf_counter())
            self.mainLogger.debug("Adding window to MDI Area, creation took %f seconds" % float(time.perf_counter() - \
                                                                                                self.localtime))
            self.mdiArea.addSubWindow(sub)

        else:
            sub.show()
            self.mdiArea.setActiveSubWindow(sub)

    # UTILITY METHODS
    def rebuildTrees(self):
        """
        At this point, the sequences and the alignments both have Window IDs applied -- but the windows
        no longer exist. Need to make sure when regenerating the windows, the old windowIDs are not lost.
        """
        for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            if node.data(role=self.SequenceRole):
                print(node.data())
                ali = {}  # empty dict needed to send to open window
                wid = node.data(role=self.WindowRole)
                seqr = node.data(role=self.SequenceRole)[0]
                ali[seqr.name] = str(seqr.seq)
                self.makeNewWindow(wid, ali, nonode=True)
                self.bioTree.setExpanded(self.bioModel.indexFromItem(node), True)

        for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
            if node.data(role=self.SequenceRole):
                print(node.data())
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
                self.projectTree.setExpanded(self.projectModel.indexFromItem(node), True)
        print(self.windows)

    def queryTrees(self):
        print("\n\nBIOROOT")
        for child in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            print("Text: ", child.text())
            print("Name: ", child.data(role=Qt.UserRole + 1))
            seqr = child.data(role=Qt.UserRole + 2)
            print("Seq: ", child.data(role=Qt.UserRole + 2))
            print("Window Index: ", child.data(role=Qt.UserRole + 3))
        print("ALIGNMENT ROOT")
        for child in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
            print("Text: ", child.text())
            print("Name: ", child.data(role=Qt.UserRole + 1))
            print("Seq: ", child.data(role=Qt.UserRole + 2))
            print("Window Index: ", child.data(role=Qt.UserRole + 3))

    def maybeClose(self):
        # TODO: CHECK IF CHANGES
        qDialog = QuitDialog(self)
        confirm = qDialog.exec()
        if confirm == 1:
            result = self.saveWorkspace()
            if result:
                print("YES")
                return True
        elif confirm == 2:
            print("NO")
            return False
        else:
            return None

    def quit(self):
        confirm = self.maybeClose()
        if confirm is not None:
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
        """ Simple method that updates the status bar process usage statistics on timer countdown"""
        mem = self.mainProcess.memory_info().rss / 1000000
        cpu = self.mainProcess.cpu_percent()
        self.memLabel.setText("CPU: " + str(cpu) + " % | RAM: " + str(round(mem, 2)) + " MB")

    def DEBUG(self):
        # STRICTLY FOR TESTING -- FAKE DATA.
        test1 = ['GPI1A', 'MSLSQDATFVELKRHVEANEKDAQLLELFEKDPARFEKFTRLFATPDGDFLFDF' +
                 'SKNRITDESFQLLMRLAKSRGVEESRNAMFSAEKINFTENRAVLHVALRNRANRP' +
                 'ILVDGKDVMPDVNRVLAHMKEFCNEIISGSWTGYTGKKITDVVNIGIGGSDLGPL' +
                 'MVTESLKNYQIGPNVHFVSNVDGTHVAEVTKKLNAETTLFIIASKTFTTQETITN' +
                 'AETAKEWFLAKAGDAGAVAKHFVALSTNVTKAVEFGIDEKNMFEFWDWVGGRYSL' +
                 'WSAIGLSIAVHIGFDNYEKLLDGAFSVDEHFVNTPLEKNIPVILAMIGVLYNNIY' +
                 'GAETHALLPYDQYMHRFAAYFQQGDMESNGKFVTRHGQRVDYSTGPIVWGEPGTN' +
                 'GQHAFYQLIHQGTRLIPADFIAPVKTLNPIRGGLHHQILLANFLAQTEALMKGKT' +
                 'AAVAEAELKSSGMSPESIAKILPHKVFEGNKPTTSIVLPVVTPFTLGALIAFYEH' +
                 'KIFVQGIIWDICSYDQWGVELGKQLAKVIQPELASADTVTSHDASTNGLIAFIKNNA']
        seq_GPI1A = Bseq.MutableSeq(test1[1], generic_protein)
        gpi1a = models.SeqR(seq_GPI1A, id=test1[0], name=test1[0])
        gpi1a.id = test1[0]
        test2 = ['GPI1B', 'MIFELFRFIFRKKKMLGYLSDLIGTLFIGDSTEKAMSLSQDATFVELKRHVEANE' +
                 'KDAQLLELFEKDPARFEKFTRLFATPDGDFLFDFSKNRITDESFQLLMRLAKSRG' +
                 'VEESRNAMFSAEKINFTENRAVLHVALRNRANRPILVDGKDVMPDVNRVLAHMKE' +
                 'FCNEIISGSWTGYTGKKITDVVNIGIGGSDLGPLMVTESLKNYQIGPNVHFVSNV' +
                 'DGTHVAEVTKKLNAETTLFIIASKTFTTQETITNAETAKEWFLAKAGDAGAVAKH' +
                 'FVALSTNVTKAVEFGIDEKNMFEFWDWVGGRYSLWSAIGLSIAVHIGFDNYEKLL' +
                 'DGAFSVDEHFVNTPLEKNIPVILAMIGVLYNNIYGAETHALLPYDQYMHRFAAYF' +
                 'QQGDMESNGKFVTRHGQRVDYSTGPIVWGEPGTNGQHAFYQLIHQGTRLIPADFI' +
                 'APVKTLNPIRGGLHHQILLANFLAQTEALMKGKTAAVAEAELKSSGMSPESIAKI' +
                 'LPHKVFEGNKPTTSIVLPVVTPFTLGALIAFYEHKIFVQGIIWDICSYDQWGVEL' +
                 'GKQLAKVIQPELASADTVTSHDASTNGLIAFIKNNA']
        seq_GPI1B = Bseq.MutableSeq(test2[1], generic_protein)
        gpi1b = models.SeqR(seq_GPI1B, id=test2[0], name=test2[0])
        test = [gpi1a, gpi1b]
        for i in test:
            self.seqInit(i)

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


class LinnaeoApp(QApplication):
    barClick = pyqtSignal()

    def __init__(self, *args):
        super().__init__(*args)
        self.installEventFilter(self)
        self.last = None
        #self.barClick.connect(self.setSizing)

    #   self._window = None

    def event(self, event):
        #print(event, event.type())
        return super().event(event)

    def eventFilter(self, obj, event):
        if event.type() in [174,175]:
            if event.type() != self.last:
                self.barClick.emit()
                self.last = event.type()
        return super().eventFilter(obj, event)

    """
    def setSizing(self):
        if self._window:
            if not self._window.beingClicked:
                print("CLICK")
                self._window.beingClicked = True
                print(self._window.beingClicked)
            elif self._window.beingClicked:
                print("UNCLICK")
                self._window.beingClicked = False
                print(self._window.beingClicked)
    """
