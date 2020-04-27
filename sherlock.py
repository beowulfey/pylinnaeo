#!/usr/bin/python3

# Bioscience components
import Bio.Seq as Bseq
from Bio.Alphabet import generic_protein
from clustalo import clustalo

# PyQt components
from PyQt5.Qt import Qt
from PyQt5.QtGui import QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import QMainWindow

# Internal components
from ui import sherlock_ui
from ui import views

# Additional libraries
import logging


def _iterTreeView(root):
    """Internal function for iterating a TreeModel"""
    def recurse(parent):
        for row in range(parent.rowCount()):
            child=parent.child(row)
            yield child
            if child.hasChildren():
                yield from recurse(child)
    if root is not None:
        yield from recurse(root)


class Sherlock(QMainWindow, sherlock_ui.Ui_MainWindow):
    """ Main Window for Sherlock App """
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)

        # Instantiation
        self.SequenceRole = Qt.UserRole + 1
        self.WindowRole = Qt.UserRole + 2
        self.alignments = {}
        self.windows = {}
        self.mainLogger = logging.getLogger("Main")
        self.bioModel = QStandardItemModel()
        self.bioRoot = QStandardItem("Sequences")
        self.projectModel = QStandardItemModel()
        self.projectRoot = QStandardItem("Project")

        # Startup functions
        self.setupUi(self)
        self.guiInit()

        self.DEBUG()

    def guiInit(self):
        """ Initialize GUI with default parameters. """
        self.bioTree.setModel(self.bioModel)
        self.bioModel.appendRow(self.bioRoot)
        self.bioTree.setExpanded(self.bioRoot.index(), True)
        self.projectTree.setModel(self.projectModel)
        self.projectModel.appendRow(self.projectRoot)
        self.projectTree.setExpanded(self.projectRoot.index(), True)

        # Slot connections
        self.bioTree.doubleClicked.connect(self.tryCreateAlignment)
        self.actionAddFolder.triggered.connect(self.tryCreateAlignment)
        self.projectTree.doubleClicked.connect(self.reopenWindow)
        #self.projectTree.dropEvent.connect(self.seqDropEvent)

    def seqDropEvent(self):
        """
        When dropping either a separate window or sequence onto either another
        sequence or another alignment window, it creates a new window using all unique components
        """

    def tryCreateAlignment(self):
        """
        This will create a new alignment from the currently selected sequences.
        Ignores any folders that were included in the selection.
        Will not duplicate alignments. Creates a window if not made.
        """
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
                self.mainLogger.debug("Detected folder ("+self.bioModel.itemFromIndex(index).text()
                                      + ") in selection --> ignoring!")
                continue

        # Check if the two sequences have been aligned before.
        # If not, align with ClustalO and create a new window from the alignment.
        seqs = list(items.values())
        if len(seqs) > 1:
            seqs.sort()
        if items and seqs not in self.alignments.values():
            wid = str(len(list(self.alignments.values())))
            self.mainLogger.debug("Storing new alignment with ID " + wid)
            self.alignments[wid] = seqs
            aligned = clustalo(items)
            self.makeNewWindow(aligned, wid)
        else:
            self.mainStatus.showMessage("Alignment already opened", msecs=5000)
            for key, value in self.alignments.items():
                if seqs == value:
                    print(key)
                    self.reopenWindow(wid=key)

    def makeNewWindow(self, ali, wid):
        sub = views.MDISubWindow()
        sub.setAttribute(Qt.WA_DeleteOnClose, False)
        widget = views.AlignSubWindow(ali)
        self.windows[wid] = sub
        sub.setWidget(widget)

        # Show window in the view panel
        self.mdiArea.addSubWindow(sub)
        node = QStandardItem(sub.windowTitle())
        node.setData(wid, self.WindowRole)
        self.projectRoot.appendRow(node)
        sub.show()

    def reopenWindow(self, wid=None):
        """
        Checks to see if a window is open already.
        If it is not, reopens the window. If it is, gives focus.
        Also refreshes the title.
        """
        print("Reopening window!--" + wid)
        if not wid:
            item = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
        else:
            for node in _iterTreeView(self.projectRoot):
                if node.data(role=self.WindowRole) == wid:
                    item = node
        try:
            sub = self.windows[item.data(role=self.WindowRole)]
            sub.setWindowTitle(item.text())
            if not sub.isVisible():
                sub.show()
            self.mdiArea.setActiveSubWindow(sub)
        except:
            pass


    # INITIAL TESTING DATA
    # Builds a basic tree model for testing.
    def DEBUG(self):
        test1 = ['GPI1A', 'MSLSQDATFVELKRHVEANEKDAQLLELFEKDPARFEKFTRLFATPDGDFLFDF'+
                 'SKNRITDESFQLLMRLAKSRGVEESRNAMFSAEKINFTENRAVLHVALRNRANRP'+
                 'ILVDGKDVMPDVNRVLAHMKEFCNEIISGSWTGYTGKKITDVVNIGIGGSDLGPL'+
                 'MVTESLKNYQIGPNVHFVSNVDGTHVAEVTKKLNAETTLFIIASKTFTTQETITN'+
                 'AETAKEWFLAKAGDAGAVAKHFVALSTNVTKAVEFGIDEKNMFEFWDWVGGRYSL'+
                 'WSAIGLSIAVHIGFDNYEKLLDGAFSVDEHFVNTPLEKNIPVILAMIGVLYNNIY'+
                 'GAETHALLPYDQYMHRFAAYFQQGDMESNGKFVTRHGQRVDYSTGPIVWGEPGTN'+
                 'GQHAFYQLIHQGTRLIPADFIAPVKTLNPIRGGLHHQILLANFLAQTEALMKGKT'+
                 'AAVAEAELKSSGMSPESIAKILPHKVFEGNKPTTSIVLPVVTPFTLGALIAFYEH'+
                 'KIFVQGIIWDICSYDQWGVELGKQLAKVIQPELASADTVTSHDASTNGLIAFIKNNA']
        seq_GPI1A = Bseq.MutableSeq(test1[1], generic_protein)
        test1alt = [test1[0], seq_GPI1A]
        test2 = ['GPI1B', 'MIFELFRFIFRKKKMLGYLSDLIGTLFIGDSTEKAMSLSQDATFVELKRHVEANE'+
                 'KDAQLLELFEKDPARFEKFTRLFATPDGDFLFDFSKNRITDESFQLLMRLAKSRG'+
                 'VEESRNAMFSAEKINFTENRAVLHVALRNRANRPILVDGKDVMPDVNRVLAHMKE'+
                 'FCNEIISGSWTGYTGKKITDVVNIGIGGSDLGPLMVTESLKNYQIGPNVHFVSNV'+
                 'DGTHVAEVTKKLNAETTLFIIASKTFTTQETITNAETAKEWFLAKAGDAGAVAKH'+
                 'FVALSTNVTKAVEFGIDEKNMFEFWDWVGGRYSLWSAIGLSIAVHIGFDNYEKLL'+
                 'DGAFSVDEHFVNTPLEKNIPVILAMIGVLYNNIYGAETHALLPYDQYMHRFAAYF'+
                 'QQGDMESNGKFVTRHGQRVDYSTGPIVWGEPGTNGQHAFYQLIHQGTRLIPADFI'+
                 'APVKTLNPIRGGLHHQILLANFLAQTEALMKGKTAAVAEAELKSSGMSPESIAKI'+
                 'LPHKVFEGNKPTTSIVLPVVTPFTLGALIAFYEHKIFVQGIIWDICSYDQWGVEL'+
                 'GKQLAKVIQPELASADTVTSHDASTNGLIAFIKNNA']
        seq_GPI1B = Bseq.MutableSeq(test2[1], generic_protein)
        test2alt = [test2[0], seq_GPI1B]
        test = [test1alt, test2alt]
        for i in list(range(0, len(test))):
            #node = models.SeqNode(test[i][0], sequence=test[i][1])
            node = QStandardItem(test[i][0])
            node.setData(test[i][1], self.SequenceRole)
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
