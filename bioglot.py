#!/usr/bin/python3

# BIOGLOT
# This is my shitty code for managing all the parts I'm going to somehow
# hack together into a piece of software, maybe.

# Bioglot components
import workspace
from ui import bioglot_ui
from ui import alignment_ui

# PyQt components
from PyQt5.QtCore import QFile, QIODevice
from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QMainWindow, QMdiSubWindow, QTextEdit

# Bioscience components
import Bio.Seq as Bseq
from Bio.Alphabet import generic_protein
import Bio.SeqUtils

# Additional libraries
import configparser
from pathlib import Path


class BioGlot(QMainWindow, bioglot_ui.Ui_MainWindow):
    """ Main Window for Bioglot App"""
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.setGui()
        self.connectSlots()
        self.DEBUG()

    def setGui(self):
        self.projectModel = QStandardItemModel()
        self.projectTree.setModel(self.projectModel)

    # GUI CONTROL: slot setup
    def connectSlots(self):
        self.projectTree.doubleClicked.connect(self.onNodeDbClick)

    # Actual slot actions go here.
    def onNodeDbClick(self, item):
        self.newSubwindow()

    # Additional UI components
    def newSubwindow(self):
        sub = QMdiSubWindow()
        # get associated sequence for click
        index = self.projectTree.selectedIndexes()[0] # can use this to select multiple!
        item = self.projectModel.itemFromIndex(index).getSeq()
        sub.setWidget(QTextEdit(str(item)))
        self.mdiArea.addSubWindow(sub)
        sub.show()

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
        root = self.projectModel.invisibleRootItem()
        for i in list(range(0, len(test))):
            node = workspace.SeqNode(test[i][0], test[i][1])
            root.appendRow(node)

        """
        # SAVED THIS FOR LATER!
        # Read in config (linux)
        config = configparser.ConfigParser()
        try:
            config.read(str(Path.home())+"/.bioglot/config.ini")
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
