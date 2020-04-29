#!/usr/bin/python3

# Bioscience components

import Bio.Seq as Bseq
from Bio.Alphabet import generic_protein
from clustalo import clustalo

# PyQt components
import PyQt5.Qt
from PyQt5.QtCore import QTimer
from PyQt5.QtGui import QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import QMainWindow, QLabel, QApplication

# Internal components
from ui import sherlock_ui
from ui import views

# Additional libraries
import sys
import os
import logging
import psutil


# TODO: Add functionality for removing sequences and alignments (from the dicts too)
# TODO: Add functionality for saving workspace.
# TODO: Add functionality for combining sequences into new alignments!


def _iterTreeView(root):
    """Internal function for iterating a TreeModel"""

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

        # Project instants
        self.windex = -1  # Acts as identifier for tracking alignments (max 2.1 billion)
        self.alignments = {}  # Alignments (stored as { windex : [name, seqs] } )
        self.windows = {}  # Windows (stored as { windex : MDISubWindow } )
        self.SequenceRole = PyQt5.Qt.Qt.UserRole + 1  # Used for storing sequence data in TreeView
        self.WindowRole = PyQt5.Qt.Qt.UserRole + 2  # Stores window ID in TreeView
        self.bioModel = QStandardItemModel()  # BioModel is shown in the top (sequence) TreeView
        self.bioRoot = QStandardItem("Folder")  # Default root node for top TreeView
        self.projectModel = QStandardItemModel()  # ProjectModel is shown in the bottom (alignment) TreeView
        self.projectRoot = QStandardItem("Folder")  # Default root node for bottom TreeView
        self.mdiArea = views.MDIArea(tabs=True)  # Create a custom MDIArea

        # Startup functions
        self.setupUi(self)  # Built by PyUic5 from my main window UI file
        self.gridLayout_2.addWidget(self.mdiArea)  # Add custom MDI area to the empty space intended to hold it
        self.guiInit()  # Additional gui setup goes here.

        self.DEBUG()

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
        self.bioTree.doubleClicked.connect(self.tryCreateAlignment)
        self.actionAlign.triggered.connect(self.tryCreateAlignment)
        self.projectTree.doubleClicked.connect(self.treeDbClick)
        # self.projectTree.dropEvent.connect(self.seqDropEvent)

    def seqDropEvent(self, event):
        """
        When dropping either a separate window or sequence onto either another
        sequence or another alignment window, it creates a new window using all unique components
        """
        print(event.source())

    def tryCreateAlignment(self, indexes=None, items=None):
        """
        This will create a new alignment from the currently selected sequences in top Tree.
        Ignores any folders that were included in the selection.
        Will not duplicate alignments. Creates a new window if alignment is new.
        "Items" input is an dictionary of {SeqName : Sequence}
        """
        if not items:
            items = {}
            for index in self.bioTree.selectedIndexes():
                # Quick and dirty way to ignore folders that are selected:
                # Only does the thing if there is a sequence present in the node.
                if self.bioModel.itemFromIndex(index).data(role=self.SequenceRole):
                    items[self.bioModel.itemFromIndex(index).text()] = \
                        str(self.bioModel.itemFromIndex(index).data(role=self.SequenceRole))
                else:
                    self.mainStatus.showMessage("Not including selected folder \"" +
                                                self.bioModel.itemFromIndex(index).text() + "\"",
                                                msecs=5000)
                    continue

        # Check if the two sequences have been aligned before.
        # If not, align with ClustalO and create a new window from the alignment.
        # If so, provide focus to that window.
        seqs = list(items.values())
        if len(seqs) > 1:
            seqs.sort()
        if items and seqs not in self.alignments.values():
            # create new unique identifier for tracking alignment!
            self.windex += 1
            wid = str(self.windex)
            self.mainLogger.debug("Alignment stored with ID " + wid + " and sequence(s) "
                                  + str([x[:10] + '..' + x[-10:] for x in seqs]))
            self.alignments[wid] = seqs

            # TODO: have this form a new thread (probably necessary for long alignments)
            # TODO: Also consider storing this as a BioPy alignment
            aligned = clustalo(items)
            self.mainLogger.debug("Sending alignment to clustal omega using default options (1 core, protein seq)")
            self.mainStatus.showMessage("Aligning selection...", msecs=1000)
            self.makeNewWindow(aligned, wid)
        else:
            self.mainStatus.showMessage("Reopening alignment!", msecs=1000)
            for key, value in self.alignments.items():
                if seqs == value:
                    self.openWindow(windowID=key)

    def makeNewWindow(self, ali, windowID):
        sub = views.MDISubWindow()
        widget = views.AlignSubWindow(ali)
        sub.setWidget(widget)

        # Show window in the view panel
        if len(ali.keys()) > 1:
            node = QStandardItem(sub.windowTitle())
            node.setData(windowID, self.WindowRole)
            self.projectRoot.appendRow(node)
        else:
            sub.setWindowTitle(list(ali.keys())[0])
        self.windows[windowID] = sub
        self.openWindow(windowID=windowID)

    def treeDbClick(self):
        print(self.projectTree.selectedIndexes()[0])
        item = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
        title = item.text()
        wid = item.data(role=self.WindowRole)
        self.openWindow(windowID=wid, title=title)

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
            if title:
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
            # node = models.SeqNode(test[i][0], sequence=test[i][1])
            node = QStandardItem(test[i][0])
            node.setData(test[i][1], self.SequenceRole)
            node.setFlags(node.flags() ^ PyQt5.Qt.Qt.ItemIsDropEnabled)
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


def main():
    logging.basicConfig(level=logging.DEBUG)  # , format="%(asctime)s:%(levelname)s:%(message)s")
    app = QApplication(sys.argv)
    try:
        with open('./ui/sherlock.sty') as f:
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
