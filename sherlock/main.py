#!/usr/bin/env python3

# Bioscience components
import copy

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
from classes import views, utilities
from ui import sherlock_ui

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
    def __init__(self, opened=None, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        # Initialize UI
        if self._instance:
            print("Closing old")
            print(self.__class___instance)
            self.__class___instance.close()
        self.setupUi(self)  # Built by PyUic5 from my main window UI file
        # System instants
        self.memLabel = QLabel()  # initiate a label for adding to status bar
        self.mainProcess = psutil.Process(os.getpid())  # process and timer are used for getting Mem/CPU usage
        self.processTimer = QTimer()  # See above. All this goes in the memLabel.
        self.mainLogger = logging.getLogger("Main")  # Logger for main window
        self.threadpool = QThreadPool()
        self.mainLogger.info("Threading with a maximum of %d threads" % self.threadpool.maxThreadCount())
        # Project instants
        self.lastClickedTree = None
        self.lastAlignment = {}
        self.titles = []  # maintains a list of sequence titles to confirm uniqueness
        self.windex = -1  # Acts as identifier for tracking alignments (max 2.1 billion)
        self.alignments = {}  # Alignments (stored as { windex : [name, seqs] } )
        self.windows = {}  # Windows (stored as { windex : MDISubWindow } )
        self.nameRole = Qt.UserRole + 1
        self.SequenceRole = Qt.UserRole + 2  # Used for storing sequence data in TreeView
        self.WindowRole = Qt.UserRole + 3  # Stores window ID in TreeView
        self.bioRoot = QStandardItem("FOLDER")
        self.bioRoot.setData("FOLDER")  # Default root node for top TreeView
        self.bioModel = views.ItemModel(  # BioModel is shown in the top (sequence) TreeView
                                        self.windows, seqTree=True)
        self.bioTree = views.TreeView()
        self.projectTree = views.TreeView()
        self.projectRoot = QStandardItem("Folder").setData("FOLDER")  # Default root node for bottom TreeView
        self.projectModel = views.ItemModel(self.windows)  # ProjectModel is shown in the bottom (alignment) TreeView
        self.mdiArea = views.MDIArea()  # Create a custom MDIArea
        self.gridLayout_2.addWidget(self.mdiArea)  # Add custom MDI area to the empty space intended to hold it

        self.guiInit()  # Additional gui setup goes here.
        self.connectSlots()
        self.DEBUG()

        if opened:
            self.bioModel = opened

    def guiInit(self):
        """ Initialize GUI with default parameters. """
        # Tree setup
        self.splitter_2.addWidget(self.bioTree)
        self.bioTree.setModel(self.bioModel)
        self.bioTree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.bioModel.appendRow(self.bioRoot)
        self.bioTree.setExpanded(self.bioModel.invisibleRootItem().index(), True)
        self.bioModel.setHorizontalHeaderLabels(["Sequences"])
        self.splitter_2.addWidget(self.projectTree)
        self.projectTree.setModel(self.projectModel)
        self.projectModel.appendRow(self.projectRoot)
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
        self.bioModel.nameCheck.connect(self.pruneNames)
        self.bioModel.nameChanged.connect(self.pruneNames)
        self.processTimer.timeout.connect(self.updateUsage)

    def instance(self, window):
        self._instance = window
    # SINGLE-USE SLOTS
    def newWorkspace(self):
        # TODO: This does not close the old window, creating a sure-fire memory leak.
        # TODO: Request save
        #self.setAttribute(Qt.WA_DeleteOnClose)
        new = Sherlock(self)
        new.instance(self)
        new.show()
        #self.hide()

    def openWorkspace(self):
        sel = QFileDialog.getOpenFileName(parent=self, caption="Open Workspace", directory=QDir.homePath(),
                                          filter="Linnaeo Workspace (*.lno);;Any (*)")
        filename = sel[0]
        print(filename)
        infile = QFile(filename)
        infile.open(QIODevice.ReadOnly)
        inF = QDataStream(infile)
        windows = {}
        newModel = views.ItemModel(windows, seqTree=True)
        self.restore_item(newModel.invisibleRootItem(),inF.readUInt32(),inF)
        #new = Sherlock(self, opened=newModel)
        infile.close()
        #new.show()

    def restore_item(self, parent, num_childs, datastream):
        try:
            print("Name: ", parent.text())
            print("DATA1: ", type(parent.data()))
            print("DATA2: ", type(parent.data(role=self.SequenceRole)))
            print("DATA3, ", type(parent.data(role=self.WindowRole)))
            print("Children: ", parent.rowCount())
        except:
            print("data not present")
        for i in range(0, num_childs):
            child = QStandardItem()
            child.read(datastream)
            num_childs = datastream.readUInt32()
            self.restore_item(child, num_childs, datastream)

    def saveWorkspace(self):
        sel = QFileDialog.getSaveFileName(parent=self, caption="Save Workspace", directory=QDir.homePath(), filter="Linnaeo Workspace (*.lno);;Any (*)")
        filename = sel[0]
        if filename[-4:] != ".lno":
            filename=str(filename+".lno")
        file = QFile(filename)
        file.open(QIODevice.WriteOnly)
        out = QDataStream(file)
        self.save_item(self.bioModel.invisibleRootItem(), out)
        file.close()

    def save_item(self, parent, datastream):
        num_childs = parent.rowCount()
        datastream.writeUInt32(num_childs)
        try:
            print("Name: ",parent.text())
            print("DATA1: ", type(parent.data()))
            print("DATA2: ", type(parent.data(role=self.SequenceRole)))
            print("DATA3, ", type(parent.data(role=self.WindowRole)))
            print("Children: ",parent.rowCount())
        except:
            print("data not present")
        for i in range(0, num_childs):
            child = parent.child(i)
            child.write(datastream)
            num_childs = child.rowCount()
            datastream.writeUInt32(num_childs)
            self.save_item(child, datastream)







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
                for line in clip:
                    if line[0] == ">" and not name:
                        print("name: ", line[:10])
                        name = line[:10]
                    elif line[0] == ">" and name:
                        print("new seq: ", line)
                        seqs.append([name, "".join(seq)])
                        self.titles.append(name)
                        name = line[:10]
                    else:
                        print("line: ", line)
                        seq.append(line.strip())
                seqs.append([name, "".join(seq)])
                for newseq in seqs:
                    self.seqCreate(newseq)
            except:
                self.mainStatus.showMessage("Please only paste in FASTA format!", msecs=1000)
            else:
                self.mainStatus.showMessage("Please choose paste destination", msecs=1000)

    def addFolder(self):
        try:
            node = self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0])
            new = QStandardItem("New Folder")
            new.setData("New Folder", role=self.nameRole)
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
                print("Nothing Selected, using last click")
                if self.lastClickedTree:
                    self.lastClickedTree.model().appendRow(new)
                else:
                    self.mainStatus.showMessage("Select a location first!", msecs=1000)

    def deleteNode(self):
        for index in self.lastClickedTree.selectedIndexes():
            node = self.lastClickedTree.model().itemFromIndex(index)
            wid = node.data(role=self.WindowRole)
            try:
                sub = self.windows[wid]
                sub.close()
                self.windows.pop(wid)
                self.alignments.pop(wid)
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
        try:
            self.mdiArea.activeSubWindow().close()
        except AttributeError:
            pass

    def closeAllTabs(self):
        self.mdiArea.closeAllSubWindows()

    # DATA AWARENESS SLOTS
    def seqDbClick(self):
        """
        This assesses what has been selected, adds that to the parent list of sequences,
        and depending on whether it is a single sequence or multiple,
        either creates, stores and (re)displays the window or first makes an alignment then does all that.
        Ignores any folders that were included in the selection.
        Will not duplicate alignments. Creates a new window only if alignment is new.
        """
        title = None
        items = {}
        # Collect the selected sequence(s)
        for index in self.bioTree.selectedIndexes():
            # Quick and dirty way to ignore folders that are selected.
            # Only does the thing if there is a sequence present in the node.
            if self.bioModel.itemFromIndex(index).data(role=self.SequenceRole):
                items[self.bioModel.itemFromIndex(index).text()] = \
                    str(self.bioModel.itemFromIndex(index).data(role=self.SequenceRole))
        seqs = list(items.values())

        # Process the sequences, and generate an alignment if >2
        # TODO: Do pairwise here if only 2!
        if len(seqs) > 1:
            # Sort the sequences to prevent duplicates and generate the alignment in a new thread.
            seqs.sort()
            worker = utilities.AlignThread(items, num_threads=self.threadpool.maxThreadCount())
            worker.start()
            worker.wait()
            aligned = worker.aligned
        else:
            aligned = items  # send single sequence

        # If this is the first time this combination has ever been selected, assign it an ID,
        # Add that window ID to the node and continue to making the window.
        if items and seqs not in self.alignments.values():
            self.windex += 1
            wid = str(self.windex)
            # TODO: Also consider storing this as a BioPy alignment
            self.alignments[wid] = seqs
            self.mainStatus.showMessage("Aligning selection...", msecs=1000) if len(seqs) > 1 else \
                self.mainStatus.showMessage("Sequence loaded", msecs=1000)
            if len(seqs) == 1:
                # Add window ID to the node
                self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0]).setData(wid, self.WindowRole)
            self.makeNewWindow(aligned, wid)
        else:
            # IF IT DOES EXIST: Messaging, and then check sequence names
            if len(seqs) > 1:
                self.mainStatus.showMessage("Reopening alignment!", msecs=1000)
            elif len(seqs) == 1:
                self.mainStatus.showMessage("Sequence loaded", msecs=1000)
                title = self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0]).text()
            else:
                # Empty else for when only a folder was selected.
                pass

            for key, value in self.alignments.items():
                # compare to previously made sequences and alignments
                if seqs == value:
                    if len(seqs) == 1:
                        print("opening single sequence")
                        self.openWindow(windowID=key, title=title)
                        if title not in self.titles:
                            self.titles.append(title)
                            self.pruneNames()
                            self.bioModel.updateNames(self.titles)
                    else:
                        self.openWindow(windowID=key)

    def alignmentDbClick(self):
        # Checks if not a folder first, then:
        # Gets the selected item (only single selection allowed), and opens the window
        item = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
        if item.data(role=self.WindowRole):
            self.openWindow(windowID=item.data(role=self.WindowRole), title=item.text())

    def dupeNameMsg(self):
        self.mainStatus.showMessage("Please choose a unique name!", msecs=1000)

    # MAIN METHODS
    def seqCreate(self, seq):
        # Input is array [name : seq]
        print(seq)
        name = seq[0]
        name = self.checkNames(name)
        node = QStandardItem(name)
        node.setData(seq[1], self.SequenceRole)
        node.setData(str(self.windex), self.WindowRole)
        self.windex += 1
        # TODO: This option allows dropping on a sequence. May consider turning off.
        node.setFlags(node.flags() ^ Qt.ItemIsDropEnabled)
        self.bioModel.appendRow(node)

    def checkNames(self, name):
        # Called upon sequence creation.
        # TODO: I NEEDED TO ADD THIS SOMEWHERE???
        # Look for ifs in TITLES!
        if name not in self.titles:
            # SAFE! You can add and return
            finalname = name
            self.titles.append(finalname)
            return finalname
        elif name[-2] == "_" and int(name[-1]):
            # if there's already a name with an _1, add a number
            # TODO: CHeck if sequence already exists!!
            finalname = str(name[:-1] + str(int(name[-1]) + 1))
            self.titles.append(finalname)
        else:
            # It's a duplicate! Better
            # TODO: CHECK SEQUENCES HERE.
            finalname = str(name + "_" + str(self.titles.count(name)))
            finalname = self.checkNames(finalname)
            self.titles.append(finalname)
        return finalname

    def makeNewWindow(self, ali, windowID):
        sub = views.MDISubWindow()
        widget = views.AlignSubWindow(ali)
        sub.setWidget(widget)
        if len(ali.keys()) > 1:
            node = QStandardItem(sub.windowTitle())
            node.setData(windowID, self.WindowRole)
            self.projectRoot.appendRow(node)
        else:
            sub.setWindowTitle(list(ali.keys())[0])
            self.titles.append(list(ali.keys())[0])
            self.bioModel.updateNames(self.titles)
        self.windows[windowID] = sub
        self.openWindow(windowID=windowID)

    def openWindow(self, windowID=None, title=None):
        """
        Checks to see if a window is open already.
        If it is not, reopens the window. If it is, gives focus.
        Also refreshes the title.
        """
        sub = self.windows[windowID]
        if title and title != sub.windowTitle():
            print("Setting new Title")
            sub.setWindowTitle(title)
        if sub:
            if sub.mdiArea() != self.mdiArea:
                print("Not in mdiArea... adding")
                self.mdiArea.addSubWindow(sub)
            else:
                sub.show()
                self.mdiArea.setActiveSubWindow(sub)

    # UTILITY METHODS
    def deselectProject(self):
        self.projectTree.clearSelection()
        self.lastClickedTree = self.bioTree

    def lastClickedSeq(self):
        self.bioModel.updateNames(self.titles)
        self.bioModel.updateLastClicked(self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0]))

    def deselectSeqs(self):
        self.bioTree.clearSelection()
        self.lastClickedTree = self.projectTree

    def pruneNames(self):
        """
        CONNECTS TO SIGNAL
        This checks which names are in "TITLES" and deletes if they were not.
        """
        names = []
        pruned = []
        for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            if node.data(self.WindowRole):
                names.append(node.text())
        for title in self.titles:
            if title not in names:
                pruned.append(title)
        print("Detected names: ", names)
        print("Stored titles: ", self.titles)
        self.titles = [x for x in self.titles and names if x not in pruned]
        print("Removed the following names: ", pruned)
        print("Now storing ", self.titles)
        self.bioModel.updateNames(self.titles)

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
        test1alt = [test1[0], seq_GPI1A]
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
        test2alt = [test2[0], seq_GPI1B]
        test = [test1alt, test2alt]
        for i in test:
            self.seqCreate(i)

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
