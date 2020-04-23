#!/usr/bin/python3

# BIOGLOT
# This is my shitty code for managing all the parts I'm going to somehow
# hack together into a piece of software, maybe.

# Bioscience components
import Bio.Seq as Bseq
from Bio.Alphabet import generic_protein
from clustalo import clustalo

# PyQt components
from PyQt5.Qt import Qt
from PyQt5.QtGui import QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import QMainWindow, QMdiSubWindow

# Bioglot components
import models
from ui import sherlock_ui
from ui import views

# Additional libraries
import logging
import configparser
from pathlib import Path

class Sherlock(QMainWindow, sherlock_ui.Ui_MainWindow):
    """ Main Window for Sherlock App"""
    def __init__(self):
        super(self.__class__, self).__init__()
        # Instantiation
        self.mainLogger = logging.getLogger("Main")
        self.bioModel = QStandardItemModel()
        self.bioRoot = QStandardItem("Sequences")
        self.projectModel = QStandardItemModel()
        self.projectRoot = QStandardItem("Project")
        # Startup functions
        self.setupUi(self)
        self.guiInit()
        self.DEBUG()

    # Initialize GUI with default parameters.
    def guiInit(self):
        # setting up trees etc.
        self.bioTree.setModel(self.bioModel)
        self.bioModel.appendRow(self.bioRoot)
        self.projectTree.setModel(self.projectModel)
        self.projectModel.appendRow(self.projectRoot)

        # slot connections go here.
        self.bioTree.doubleClicked.connect(self.newSubWindow)
        self.projectTree.doubleClicked.connect(self.reopenWindow)
        self.actionAddFolder.triggered.connect(self.newSubWindow)

    # Slot actions go here
    def newSubWindow(self):
        """
        This will create a new alignment window from whatever the selected sequences were.
        Ignores any folders that were included in the selection.
        """
        # get associated sequence for click
        indexes = self.bioTree.selectedIndexes()
        items = {}
        for index in indexes:
            # Quick and dirty way to ignore folders that are selected too.
            try:
                items[self.bioModel.itemFromIndex(index).text()] = \
                    str(self.bioModel.itemFromIndex(index).getSeq())
            except AttributeError:
                self.mainStatus.showMessage("Not including selected folder \"" +
                                            self.bioModel.itemFromIndex(index).text() + "\"",
                                            msecs=5000)
                self.mainLogger.debug("Detected folder ("+self.bioModel.itemFromIndex(index).text()
                                      + ") in selection --> ignoring!")
                continue
        # Align with ClustalO and create a new window from the alignment.
        aligned = clustalo(items)
        sub = views.MDISubWindow()
        sub.setAttribute(Qt.WA_DeleteOnClose, False)

        sub.setWidget(views.AlignSubWindow(aligned))
        self.mdiArea.addSubWindow(sub)
        node = models.WorkspaceNode("New Window", sub)
        self.projectRoot.appendRow(node)
        sub.show()

    def reopenWindow(self):
        """
        Checks to see if a window is open already.
        If it is not, reopens the window. If it is, gives focus.
        """
        item = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
        try:
            sub = item.getWindow()
            sub.setWindowTitle(item.text())
            self.mdiArea.setActiveSubWindow(sub)
            if not self.mdiArea.isVisible(sub):
                sub.show()
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
        seq_GPI1A = Bseq.MutableSeq(test1[1],generic_protein)
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
        seq_GPI1B = Bseq.MutableSeq(test2[1],generic_protein)
        test2alt = [test2[0], seq_GPI1B]
        test = [test1alt, test2alt]
        for i in list(range(0, len(test))):
            node = models.SeqNode(test[i][0], test[i][1])
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
