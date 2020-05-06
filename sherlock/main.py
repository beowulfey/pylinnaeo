#!/usr/bin/env python3

# Bioscience components
import copy

import Bio
import Bio.Seq as Bseq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as Bseqio
from Bio.Alphabet import generic_protein
from clustalo import clustalo

# PyQt components
from PyQt5.Qt import Qt
from PyQt5.QtCore import QThreadPool, QTimer, QDir, QFile, QIODevice, QDataStream
from PyQt5.QtGui import QStandardItem
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QAbstractItemView, QShortcut, QFileDialog

# Internal components
from sherlock.classes import models, views, utilities
from sherlock.ui import sherlock_ui

import sys
import os
import logging
import psutil


# TODO: #1 Fix issues with renaming sequences: trying to make it so they cannot be duplicated
# TODO: Add functionality for removing sequences and alignments (from the dicts too)
# TODO: Add functionality for saving workspace.
# TODO: Add functionality for combining sequences into new alignments!




class Sherlock(QMainWindow, sherlock_ui.Ui_MainWindow):
    """
    Constructor for the Main Window of the Sherlock App
    Contains all the user interface functions, as well as underlying code for major features.
    The bulk of the program is located here.
    """

    _instance = None

    def __init__(self, *args, trees=None, data=None):

        super(self.__class__, self).__init__(*args)
        # Initialize UI
        self.setupUi(self)  # Built by PyUic5 from my main window UI file

        # System instants; status bar and threads
        self.memLabel = QLabel()
        self.mainProcess = psutil.Process(os.getpid())
        self.processTimer = QTimer()
        self.mainLogger = logging.getLogger("Main")
        self.threadpool = QThreadPool()
        self.mainLogger.info("Threading with a maximum of %d threads" % self.threadpool.maxThreadCount())

        # Project instants; inherent variables for logic.
        self.lastClickedTree = None
        self.lastAlignment = {}
        self.SequenceRole = Qt.UserRole + 2
        self.WindowRole = Qt.UserRole + 3

        # MDI Window
        self.mdiArea = views.MDIArea()
        self.gridLayout_2.addWidget(self.mdiArea)

        # Tree stuff
        self.bioTree = views.TreeView()
        self.projectTree = views.TreeView()
        self.windows = {}  # Windows stored as { windex : MDISubWindow }
        self.windex = 0  # Acts as identifier for tracking alignments (max 2.1 billion)

        if trees:
            # For loading a window on File>Open
            print("Using saved workspace")
            self.bioModel, self.projectModel = trees
            self.sequences, self.windex, self.alignments = data
            self.bioRoot = self.bioModel.invisibleRootItem().child(0)
            self.projectRoot = self.projectModel.invisibleRootItem().child(0)
            self.windows = {}
            self.sequences = {}
            self.titles = self.sequences.keys()
            self.rebuildTrees()
            print("TITLES: ", str(self.titles))
            print("ALIGNS: ", str(self.alignments))
            print("WINDEX: ", str(self.windex))

        else:
            self.sequences = {}  # Stored as {WindowID:SeqRecord (or [SeqRecords]) }
            self.alignments = {}  # Alignments (stored as { windex : [name, seqs] } )
            self.titles = []  # maintains a list of sequence titles to confirm uniqueness
            self.bioRoot = QStandardItem("Folder")
            self.bioModel = views.ItemModel(self.windows, seqTree=True)
            self.bioRoot.setData("Folder")
            self.bioModel.appendRow(self.bioRoot)
            self.DEBUG()
            self.projectRoot = QStandardItem("Folder")
            self.projectModel = views.ItemModel(self.windows)
            self.projectModel.appendRow(self.projectRoot)

        # Other functions
        self.guiInit()
        self.connectSlots()

    # TODO: I need to make it so SEQUENCES ARE SAVED
    # TODO: I need to then make it so I REBUILD THE ALIGNMENT using that data
    # That way, I can save exactly which alignment has which baseline s equences.
    # May require some SIGNIFICANT RESTRUCTURING.

    def rebuildTrees(self):
        seqs = []
        for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
            ali = node.data(role=self.SequenceRole)
            #self.windex = self.windex + 1
            #wid = self.windex
            #node.setData(wid, role=self.WindowRole)
            for key, value in self.alignments.items():
                if value == ali:
                    node.setData(key, role=self.WindowRole)
                    self.windex = int(node.data(role=self.WindowRole))+1
                    print("NEW WINDEX: ", self.windex)

                worker = utilities.AlignThread(seqs, num_threads=self.threadpool.maxThreadCount())
                worker.start()
                worker.wait()
                ali = worker.aligned
                self.makeNewWindow(ali, node.data(role=self.WindowRole), node)
                self.mdiArea.closeAllSubWindows()
        #for node in nodes:
        #    self.projectModel.removeRow(node.index().row(), node.index().parent())

    def queryTrees(self):
        print("BIOROOT")
        for child in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            print("Text: ", child.text())
            print("Name: ", child.data(role=Qt.UserRole + 1))
            print("Seq: ", child.data(role=Qt.UserRole + 2))
            print("Window Index: ", child.data(role=Qt.UserRole + 3))
        print("ALIGNMENT ROOT")
        for child in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
            print("Text: ", child.text())
            print("Name: ", child.data(role=Qt.UserRole + 1))
            print("Seq: ", child.data(role=Qt.UserRole + 2))
            print("Window Index: ", child.data(role=Qt.UserRole + 3))

    def guiInit(self):
        """ Initialize GUI with default parameters. """
        # Tree setup
        self.splitter_2.addWidget(self.bioTree)
        self.bioTree.setModel(self.bioModel)
        self.bioTree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.bioTree.setExpanded(self.bioModel.invisibleRootItem().index(), True)
        self.bioModel.setHorizontalHeaderLabels(["Sequences"])
        self.splitter_2.addWidget(self.projectTree)
        self.projectTree.setModel(self.projectModel)
        self.projectTree.setExpanded(self.projectModel.invisibleRootItem().index(), True)
        self.projectModel.setHorizontalHeaderLabels(["Alignments"])

        # Status bar setup
        self.updateUsage()
        self.statusBar().addPermanentWidget(self.memLabel)
        self.processTimer.setInterval(1000)
        self.processTimer.start()

    def connectSlots(self):
        # Toolbar and MenuBar
        # FILE
        self.actionNew.triggered.connect(self.newWorkspace)
        self.actionOpen.triggered.connect(self.openWorkspace)
        self.actionSave.triggered.connect(self.saveWorkspace)

        # EDIT
        self.actionCopy.triggered.connect(self.copyOut)
        self.actionPaste.triggered.connect(self.pasteInto)
        #self.actionPreferences.triggered.connect(self.openPrefWindow)

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

    # SINGLE-USE SLOTS
    def newWorkspace(self):
        # TODO: This does not close the old window, creating a sure-fire memory leak.
        # TODO: Request save
        self.setAttribute(Qt.WA_DeleteOnClose)
        new = Sherlock(self)
        new.instance(self)
        new.show()
        self.hide()

    def openWorkspace(self):
        sel = QFileDialog.getOpenFileName(parent=self, caption="Open Workspace", directory=QDir.homePath(),
                                          filter="Linnaeo Workspace (*.lno);;Any (*)")
        filename = sel[0]
        self.mainLogger.debug("Opening file: "+str(filename))
        file = QFile(filename)
        file.open(QIODevice.ReadOnly)
        fstream = QDataStream(file)
        titles = fstream.readQVariantList()
        windex = fstream.readUInt32()
        windows = {}
        aligns = fstream.readQVariantHash()
        newModel = views.ItemModel(windows, seqTree=True)
        newProjectModel = views.ItemModel(windows)
        self.restore_item(newModel.invisibleRootItem(), fstream)
        self.restore_item(newProjectModel.invisibleRootItem(), fstream)
        new = Sherlock(self, trees=[newModel, newProjectModel], data=[titles, windex, aligns])
        file.close()
        new.show()
        self.hide()

    def restore_item(self, parent, datastream, num_childs=None):
        print("Starting restore")
        if not num_childs:
            print("First line: reading UInt32")
            num_childs = datastream.readUInt32()
            print(num_childs)
        for i in range(0, num_childs):
            print("Reading node")
            child = QStandardItem()
            child.read(datastream)
            parent.appendRow(child)
            print(child.data(role=Qt.UserRole+1))
            print(child.data(role=Qt.UserRole+2))
            print(child.data(role=Qt.UserRole+3))
            num_childs = datastream.readUInt32()
            if num_childs > 0:
                print("reading children")
                self.restore_item(child, datastream, num_childs)

    def saveWorkspace(self):
        sel = QFileDialog.getSaveFileName(parent=self, caption="Save Workspace", directory=QDir.homePath(),
                                          filter="Linnaeo Workspace (*.lno);;Any (*)")
        filename = sel[0]
        if filename[-4:] != ".lno":
            filename = str(filename+".lno")
        file = QFile(filename)
        file.open(QIODevice.WriteOnly)
        out = QDataStream(file)

        self.mainLogger.debug("Beginning file save")
        out.writeQVariantList(self.titles)
        out.writeUInt32(self.windex)
        out.writeQVariantHash(self.alignments)
        #####################################
        # Sequence View
        ###
        out.writeUInt32(self.bioModel.invisibleRootItem().rowCount())
        for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            node.write(out)
            out.writeUInt32(node.rowCount())
        out.writeUInt32(self.projectModel.invisibleRootItem().rowCount())
        for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
            print("Text: ", str(node.text()))
            print("Data1: ", str(node.data()))
            print("Data2: ", str(node.data(role=self.SequenceRole)))
            print("Data3: ", str(node.data(role=self.WindowRole)))
            node.write(out)
            out.writeUInt32(node.rowCount())
        self.mainLogger.debug("Save complete")
        #####################################
        # Alignment View
        # Does not save any metadata! Only the two sequences
        file.flush()
        file.close()

    def copyOut(self):
        if self.lastClickedTree == self.bioTree:
            seqs = []
            indices = self.lastClickedTree.selectedIndexes()
            for index in indices:
                node = self.bioModel.itemFromIndex(index)
                seqs.append(SeqRecord(node.data(role=self.SequenceRole), id=node.text()).format("fasta"))
            QApplication.clipboard().setText("".join(seqs))
            if len(seqs) > 1:
                self.mainStatus.showMessage("Copied sequences to clipobard!", msecs=1000)
            elif len(seqs) == 1:
                self.mainStatus.showMessage("Copied sequence to clipobard!", msecs=1000)

        elif self.lastClickedTree == self.projectTree:
            self.mainStatus.showMessage("Please use File>Export to save alignments!", msecs=1000)

    def pasteInto(self):
        # FASTA DETECTION
        if self.lastClickedTree == self.bioTree:
            name = None
            seq = []
            seqs = []
            try:
                clip = str(QApplication.clipboard().text()).splitlines()
                if clip[0][0] == ">":
                    for line in clip:
                        if line[0] == ">" and not name:
                            #self.mainLogger.debug(str(line[:10]))
                            self.pruneNames()
                            name = line[:10]

                        elif line[0] == ">" and name:

                            #self.mainLogger.debug(str(line))
                            seqs.append([name, "".join(seq)])
                            self.pruneNames()
                            name = line[:10]
                        else:
                            #self.mainLogger.debug(str(line))
                            seq.append(line.strip())
                    # TODO: Convert this into SequenceRecord!!
                    seqs.append([name, Bseq.MutableSeq("".join(seq))])
                    for newseq in seqs:
                        self.seqInit(newseq)
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
                title = sub.windowTitle()
                self.mainLogger.debug("Deleting node from tree: "+str(sub.windowTitle()))
                sub.close()
                print(self.titles)
                self.titles.remove(title)
                print(self.titles)
                self.windows.pop(wid)
                self.alignments.pop(wid)
                self.pruneNames()
                self.bioModel.updateNames(self.titles)
                self.bioModel.updateWindows(self.windows)
                self.projectModel.updateNames(self.titles)
                self.projectModel.updateWindows(self.windows)
            except KeyError:
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
        seqs = list(seqarray.values())
        if len(seqs) > 1:
            # Sort the sequences to prevent duplicates and generate the alignment in a new thread.
            seqs.sort()
            worker = utilities.AlignThread(seqarray, num_threads=self.threadpool.maxThreadCount())
            worker.start()
            worker.wait()
            aligned = worker.aligned
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
        items = {}
        combo = []
        # Collect the selected sequence(s)
        for index in self.bioTree.selectedIndexes():
            # Need to make a dictionary of { Name:Sequence } for sending to ClustalO.
            # Only does the thing if there is a sequence present in the node.
            if self.bioModel.itemFromIndex(index).data(role=self.SequenceRole):
                seqr = self.bioModel.itemFromIndex(index).data(role=self.SequenceRole)[0]
                items[seqr.id] = str(seqr.seq)
                combo.append(seqr)
        combo.sort()
        aligned = self.callAlign(items)

        if items:
            if combo in self.sequences.values():
                for key, value in self.sequences.items():
                    print("COMBO: ",combo)
                    print("VALUE: ", value)
                    if combo == value:
                        wid = key
                        try:
                            sub = self.windows[wid]
                            print("Reopening saved window")
                            self.openWindow(sub)
                        except KeyError:
                            print("Opening new window")
                            self.makeNewWindow(wid, aligned)
            else:
                wid = self.windex + 1
                self.sequences[wid] = combo
                self.makeNewWindow(wid, aligned)

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
        # Input is SeqRecord. Give it a unique WID and add to list of all SEQs.
        wid = self.windex + 1
        sid = seqr.id
        # Check if title already exists, and if so, changes it.
        sid, self.titles = utilities.checkName(sid, self.titles)
        if sid != seqr.id:
            seqr.id = sid
        print("SEQ INIT: ", sid)
        print("Updated titles: ", self.titles)
        # Adds to the list of sequences, by its Window ID
        self.sequences[wid] = [seqr]
        print("Updated sequences: \n", self.sequences)
        node = QStandardItem(sid)
        node.setData(sid)
        node.setData([seqr], self.SequenceRole)
        node.setData(wid, self.WindowRole)
        node.setFlags(node.flags() ^ Qt.ItemIsDropEnabled)
        self.bioModel.appendRow(node)
        self.windex = wid

    def makeNewWindow(self, wid, ali):
        print(ali)
        sub = views.MDISubWindow()
        widget = views.AlignSubWindow(ali)
        sub.setWidget(widget)
        print(sub)
        if len(ali.keys()) > 1:
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
            sub.setWindowTitle(seq.id)
            self.windows[wid] = sub
            self.bioModel.updateNames(self.titles)
            self.bioModel.updateWindows(self.windows)
        self.openWindow(sub)

    def openWindow(self, sub):
        """
        Checks to see if a window is open already.
        If it is not, reopens the window. If it is, gives focus.
        Also refreshes the title.
        """
        self.queryTrees()
        if sub.mdiArea() != self.mdiArea:
            self.mainLogger.debug("Adding window to MDI Area")
            self.mdiArea.addSubWindow(sub)
        else:
            sub.show()
            self.mdiArea.setActiveSubWindow(sub)

    # UTILITY METHODS
    def instance(self, window):
        """Saves instance upon loading a new window"""
        # TODO: DEBUG THIS PROCESS!
        self._instance = window

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
        self.pruneNames()
        self.bioModel.updateNames(self.titles)

    def pruneNames(self):
        """
        CONNECTS TO SIGNAL
        This checks which names are in "TITLES" and deletes if they were not.
        """
        names = []
        pruned = []
        for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            if node.data(self.SequenceRole):
                names.append(node.text())
        for title in self.titles:
            if title not in names:
                pruned.append(title)
        self.mainLogger.debug("Tree names: "+str(names)+" vs. Stored names: "+str(self.titles))
        self.titles = [x for x in self.titles and names if x not in pruned]
        self.mainLogger.debug("Removed "+str(pruned)+", leaving "+str(self.titles))

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
        gpi1a = models.SeqR(seq_GPI1A, id=test1[0])
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
        gpi1b = models.SeqR(seq_GPI1B, id=test2[0])
        test = [gpi1a, gpi1b]
        for i in test:
            self.seqInit(i)

        """
        # SAVED THIS FOR LATER!
        # Read in config (linux)
        config = configparser.ConfigParser()
        try:
            config.read(str(Path.home())+"/.sherlock/config.ini")
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


def main():
    logging.basicConfig(level=logging.DEBUG)  # , format="%(asctime)s:%(levelname)s:%(message)s")
    app = QApplication(sys.argv)

    try:
        with open('ui/sherlock.sty') as f:
            print("Read in stylesheet")
            style = f.read()
    except IOError:
        print('Could not read stylesheet.')
        style = ""

    if style:
        app.setStyleSheet(style)
    window = Sherlock()
    window.show()

    app.exec_()


if __name__ == '__main__':
    main()
