#!/usr/bin/env python3

# Bioscience components
import Bio.Seq as Bseq
from Bio.Alphabet import generic_protein
from clustalo import clustalo

# PyQt components
from PyQt5.QtCore import QThreadPool, QTimer

from PyQt5.Qt import Qt
from PyQt5.QtGui import QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel

# Internal components
from sherlock.classes import views, sherlock_ui

import sys
import os
import logging
import psutil


# TODO: #1 Fix issues with renaming sequences: trying to make it so they cannot be duplicated
# TODO: Add functionality for removing sequences and alignments (from the dicts too)
# TODO: Add functionality for saving workspace.
# TODO: Add functionality for combining sequences into new alignments!


def iterTreeView(root):
    """
    Internal function for iterating a TreeModel.
    Usage: for node in _iterTreeView(root): etc.
    """
    def recurse(parent):
        for row in range(parent.rowCount()):
            child = parent.child(row)
            yield child
            if child.hasChildren():
                yield from recurse(child)
    if root is not None:
        yield from recurse(root)


class Sherlock(QMainWindow, sherlock_ui.Ui_MainWindow):
    """
    Constructor for the Main Window of the Sherlock App
    Contains all the user interface functions, as well as underlying code for major features.
    The bulk of the program is located here.
    """
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        # System instants
        self.memLabel = QLabel()  # initiate a label for adding to status bar
        self.mainProcess = psutil.Process(os.getpid())  # process and timer are used for getting Mem/CPU usage
        self.processTimer = QTimer()  # See above. All this goes in the memLabel.
        self.mainLogger = logging.getLogger("Main")  # Logger for main window

        # Startup functions
        self.setupUi(self)  # Built by PyUic5 from my main window UI file

        # Project instants
        self.titles = []  # maintains a list of sequence titles to confirm uniqueness
        self.windex = -1  # Acts as identifier for tracking alignments (max 2.1 billion)
        self.alignments = {}  # Alignments (stored as { windex : [name, seqs] } )
        self.windows = {}  # Windows (stored as { windex : MDISubWindow } )
        self.SequenceRole = Qt.UserRole + 1  # Used for storing sequence data in TreeView
        self.WindowRole = Qt.UserRole + 2  # Stores window ID in TreeView
        self.bioRoot = QStandardItem("Folder")  # Default root node for top TreeView
        self.bioModel = views.ItemModel(self.bioRoot,  # BioModel is shown in the top (sequence) TreeView
                                        self.windows, seqs=self.alignments, names=self.titles)
        self.projectRoot = QStandardItem("Folder")  # Default root node for bottom TreeView
        self.projectModel = views.ItemModel(self.projectRoot, self.windows)  # ProjectModel is shown in the bottom (alignment) TreeView
        self.mdiArea = views.MDIArea(tabs=True)  # Create a custom MDIArea
        self.gridLayout_2.addWidget(self.mdiArea)  # Add custom MDI area to the empty space intended to hold it

        self.threadpool = QThreadPool()
        self.mainLogger.debug("Multithreading with maximum %d threads" % self.threadpool.maxThreadCount())

        self.guiInit()  # Additional gui setup goes here.

    def guiInit(self):
        """ Initialize GUI with default parameters. """
        # Tree setup
        self.bioTree.setModel(self.bioModel)
        self.bioModel.appendRow(self.bioRoot)
        self.bioTree.setExpanded(self.bioRoot.index(), True)
        self.bioModel.setHorizontalHeaderLabels(["Sequences"])
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
        self.actionAlign.triggered.connect(self.seqDbClick)
        self.projectTree.doubleClicked.connect(self.alignmentDbClick)

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
        for index in self.bioTree.selectedIndexes():
            # Quick and dirty way to ignore folders that are selected.
            # Only does the thing if there is a sequence present in the node.
            if self.bioModel.itemFromIndex(index).data(role=self.SequenceRole):
               items[self.bioModel.itemFromIndex(index).text()] = \
                    str(self.bioModel.itemFromIndex(index).data(role=self.SequenceRole))
            else:
                self.mainStatus.showMessage("Not including selected folder \"" +
                                            self.bioModel.itemFromIndex(index).text() + "\"",
                                            msecs=1000)
        # Check if the sequence(s) have been aligned before.
        # If not, align with ClustalO and create a new window from the alignment.
        # If so, provide focus to that window.
        seqs = list(items.values())
        if len(seqs) > 1:
            seqs.sort()

            aligned = clustalo(items)
        if items and (seqs not in self.alignments.values()):
            # create new unique identifier for tracking alignment or sequence!
            self.windex += 1
            wid = str(self.windex)
            self.mainLogger.debug("Alignment stored with ID " + wid + " and sequence(s) "
                                  + str([x[:10] + '..' + x[-10:] for x in seqs]))

            self.alignments[wid] = seqs

            # TODO: have this form a new thread (probably necessary for long alignments)
            # TODO: Also consider storing this as a BioPy alignment

            self.mainLogger.debug("Sending alignment to clustal omega using default options (1 core, protein seq)")
            self.mainStatus.showMessage("Aligning selection...", msecs=1000) if len(seqs) > 1 else \
                self.mainStatus.showMessage("Sequence loaded", msecs=1000)
            self.makeNewWindow(aligned, wid)
        else:
            if len(seqs) > 1:
                self.mainStatus.showMessage("Reopening alignment!", msecs=1000)
            else:
                self.mainStatus.showMessage("Sequence loaded", msecs=1000)
                title = self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0]).text()
                print(title)
            for key, value in self.alignments.items():
                if seqs == value:
                    if len(seqs) == 1:
                        # TODO #1: THIS CURRENTLY RENAMES BECAUSE IT'S NOT DETECTING IF A NAME WAS UNCHANGED!
                        # REDO THIS WHOLE BLOCK
                        if title and title in self.titles:
                            print("duplicate title")
                            title = title+"_"+str(self.titles.count(title))
                            self.titles.append(title)
                        self.openWindow(windowID=key, title=title)
                        for node in iterTreeView(self.bioRoot):
                            if node.data(role=self.WindowRole) == key:
                                node.setText(title)
                        self.pruneNames()
                        # TO HERE!
                    else:
                        self.openWindow(windowID=key)
                        self.pruneNames()

    def alignmentDbClick(self):
        # Checks if not a folder first, then:
        # Gets the selected item (only single selection allowed), and opens the window
        try:
            item = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
            self.openWindow(windowID=item.data(role=self.WindowRole), title=item.text())
        except:     # TODO: Make specific
            print("Not an alignment")

    def pruneNames(self):
        # TODO: Change to updateTrees and generalize
        # This checks which names are in there and deletes if they were not.
        names = []
        pruned = []
        for node in iterTreeView(self.bioRoot):
            names.append(node.text())
        for title in self.titles:
            if title not in names:
                pruned.append(title)
        print("Detected names: ", names)
        print("Stored titles: ", self.titles)
        self.titles = [x for x in self.titles and names if x not in pruned]
        print("Removed the following names: ", pruned)
        print("Now storing ",  self.titles)

    def checkDropType(self):
        pass

    def alignCreate(self):
        pass
# Input is array of sequence arrays, each sub array is [name : seq]


    def seqCreate(self, sequence):
        pass


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
        self.windows[windowID] = sub
        self.openWindow(windowID=windowID)

    def openWindow(self, windowID=None, title=None):
        """
        Checks to see if a window is open already.
        If it is not, reopens the window. If it is, gives focus.
        Also refreshes the title.
        """
        sub = None
        if isinstance(windowID, str):
            # if WindowID is a string, that means it was sent by a deliberate search.
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
        for i in list(range(0, len(test))):
            node = QStandardItem(test[i][0])
            node.setData(test[i][1], self.SequenceRole)
            node.setData(str(self.windex), self.WindowRole)
            self.windex += 1
            node.setFlags(node.flags() ^ Qt.ItemIsDropEnabled)
            self.bioRoot.appendRow(node)

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

