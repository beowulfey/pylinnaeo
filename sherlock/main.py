#!/usr/bin/env python3

# Bioscience components
import Bio.Seq as Bseq
from Bio.Alphabet import generic_protein
from clustalo import clustalo

# PyQt components
from PyQt5.Qt import Qt
from PyQt5.QtCore import QThreadPool, QTimer
from PyQt5.QtGui import QStandardItem
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QAbstractItemView

# Internal components
from sherlock.classes import views, utilities
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

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        # Initialize UI
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
        self.SequenceRole = Qt.UserRole + 1  # Used for storing sequence data in TreeView
        self.WindowRole = Qt.UserRole + 2  # Stores window ID in TreeView
        self.bioRoot = QStandardItem("Folder")  # Default root node for top TreeView
        self.bioModel = views.ItemModel(self.bioRoot,  # BioModel is shown in the top (sequence) TreeView
                                        self.windows, seqTree=True)
        self.bioTree = views.TreeView()
        self.projectTree = views.TreeView()
        self.projectRoot = QStandardItem("Folder")  # Default root node for bottom TreeView
        self.projectModel = views.ItemModel(self.projectRoot,
                                            self.windows)  # ProjectModel is shown in the bottom (alignment) TreeView
        self.mdiArea = views.MDIArea()  # Create a custom MDIArea
        self.gridLayout_2.addWidget(self.mdiArea)  # Add custom MDI area to the empty space intended to hold it

        self.guiInit()  # Additional gui setup goes here.
        self.DEBUG()

    def guiInit(self):
        """ Initialize GUI with default parameters. """
        # Tree setup
        self.splitter_2.addWidget(self.bioTree)
        self.bioTree.setModel(self.bioModel)
        self.bioTree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.bioModel.appendRow(self.bioRoot)
        self.bioTree.setExpanded(self.bioRoot.index(), True)
        self.bioModel.setHorizontalHeaderLabels(["Sequences"])
        self.splitter_2.addWidget(self.projectTree)
        self.projectTree.setModel(self.projectModel)
        self.projectModel.appendRow(self.projectRoot)
        self.projectTree.setExpanded(self.projectRoot.index(), True)
        self.projectModel.setHorizontalHeaderLabels(["Alignments"])

        # Status bar setup
        self.updateUsage()
        self.statusBar().addPermanentWidget(self.memLabel)
        self.processTimer.setInterval(1000)
        self.processTimer.start()

        # Slot connections
        self.processTimer.timeout.connect(self.updateUsage)
        self.bioTree.doubleClicked.connect(self.seqDbClick)
        self.bioTree.generalClick.connect(self.deselectProject)
        self.bioModel.nameCheck.connect(self.pruneNames)
        self.bioModel.nameChanged.connect(self.pruneNames)
        self.bioModel.dupeName.connect(self.dupeNameMsg)
        self.projectTree.doubleClicked.connect(self.alignmentDbClick)
        self.projectTree.generalClick.connect(self.deselectSeqs)

        # Toolbar and MenuBar
        self.actionAlign.triggered.connect(self.seqDbClick)
        self.actionNewFolder.triggered.connect(self.addFolder)
        self.actionTile.triggered.connect(self.tileWindows)
        self.actionCascade.triggered.connect(self.cascadeWindows)
        # TODO: Set this is an editable preference
        self.actionToggle_Tabs.triggered.connect(self.mdiArea.toggleTabs)
        self.actionClose_all.triggered.connect(self.closeTabs)
        self.actionDelete.triggered.connect(self.deleteNode)

    # SINGLE-USE SLOTS
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

    def closeTabs(self):
        self.mdiArea.closeAllSubWindows()

    def deselectProject(self):
        self.projectTree.clearSelection()
        self.lastClickedTree = self.bioTree
        print(self.lastClickedTree.selectedIndexes())

    def deselectSeqs(self):
        self.bioTree.clearSelection()
        self.lastClickedTree = self.projectTree

    def addFolder(self):
        try:
            node = self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0])
            if not node.data(role=self.WindowRole):
                node.appendRow(QStandardItem("New Folder"))
            else:
                self.bioModel.appendRow(QStandardItem("New Folder"))
        except IndexError:
            try:
                node = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
                if not node.data(role=self.WindowRole):
                    node.appendRow(QStandardItem("New Folder"))
                else:
                    self.projectModel.appendRow(QStandardItem("New Folder"))
            except IndexError:
                print("Nothing Selected, using last click")
                if self.lastClickedTree:
                    self.lastClickedTree.model().appendRow(QStandardItem("New Folder"))
                else:
                    self.mainStatus.showMessage("Select a location first!", msecs=1000)

    def tileWindows(self):
        self.mdiArea.tileSubWindows()

    def cascadeWindows(self):
        self.mdiArea.cascadeSubWindows()

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
        try:
            item = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
            self.openWindow(windowID=item.data(role=self.WindowRole), title=item.text())
        except:  # TODO: Make specific
            print("Not an alignment")

    def dupeNameMsg(self):
        self.mainStatus.showMessage("Please choose a unique name!", msecs=1000)

    # MAIN METHODS
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
                self.mdiArea.addSubWindow(sub)
            else:
                sub.show()
                self.mdiArea.setActiveSubWindow(sub)

    def checkDropType(self):
        pass

    def seqCreate(self, seqs):
        # Input is array of sequence arrays, each sub array is [name : seq]
        for i in list(range(0, len(seqs))):
            node = QStandardItem(seqs[i][0])
            node.setData(seqs[i][1], self.SequenceRole)
            node.setData(str(self.windex), self.WindowRole)
            self.windex += 1
            # TODO: This option allows dropping on a sequence. May consider turning off.
            node.setFlags(node.flags() ^ Qt.ItemIsDropEnabled)
            self.bioModel.appendRow(node)

    # UTILITY METHODS
    def pruneNames(self):
        # This checks which names are in there and deletes if they were not.
        names = []
        pruned = []
        for node in utilities.iterTreeView(self.bioRoot):
            names.append(node.text())
        for title in self.titles:
            if title not in names:
                pruned.append(title)
        # print("Detected names: ", names)
        # print("Stored titles: ", self.titles)
        self.titles = [x for x in self.titles and names if x not in pruned]
        # print("Removed the following names: ", pruned)
        # print("Now storing ", self.titles)
        self.bioModel.updateNames(self.titles)

    def updateUsage(self):
        """ Simple method that updates the status bar process usage statistics on timer countdown"""
        mem = self.mainProcess.memory_info().rss / 1000000
        cpu = self.mainProcess.cpu_percent()
        self.memLabel.setText("CPU: " + str(cpu) + " % | RAM: " + str(round(mem, 2)) + " MB")

    def DEBUG(self):
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
        self.seqCreate(test)

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
        with open('sherlock.sty') as f:
            print("Read in stylesheet")
            style = f.read()
    except IOError:
        print('Could not read stylesheet.')
        style = ""

    # if style:
    #    app.setStyleSheet(style)
    window = Sherlock()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
