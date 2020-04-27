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
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        # Instantiation
        self.alignments = []
        self.windows = {}
        self.mainLogger = logging.getLogger("Main")
        self.bioModel = QStandardItemModel()
        self.bioModel.setItemPrototype(models.SeqNode())
        self.bioRoot = models.SeqNode("Sequences")
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
        self.bioTree.setExpanded(self.bioRoot.index(), True)
        self.projectTree.setModel(self.projectModel)
        self.projectModel.appendRow(self.projectRoot)
        self.projectTree.setExpanded(self.projectRoot.index(), True)

        # slot connections
        self.bioTree.doubleClicked.connect(self.tryCreateAlignment)
        self.projectTree.doubleClicked.connect(self.reopenWindow)
        self.actionAddFolder.triggered.connect(self.tryCreateAlignment)

    def tryCreateAlignment(self):
        """
        This will create a new alignment from whatever the selected sequences were.
        Ignores any folders that were included in the selection.
        Will not duplicate alignments.
        """
        items = {}
        for index in self.bioTree.selectedIndexes():
            # Quick and dirty way to ignore folders that are selected too.
            if self.bioModel.itemFromIndex(index).sequence():
                items[self.bioModel.itemFromIndex(index).text()] = \
                    str(self.bioModel.itemFromIndex(index).sequence())
            else:
                self.mainStatus.showMessage("Not including selected folder \"" +
                                            self.bioModel.itemFromIndex(index).text() + "\"",
                                            msecs=5000)
                self.mainLogger.debug("Detected folder ("+self.bioModel.itemFromIndex(index).text()
                                      + ") in selection --> ignoring!")
                continue

        # Check if the two sequences have been aligned before.
        # If not, align with ClustalO and create a new window from the alignment.
        if items and list(items.values()) not in self.alignments:
            aligned = clustalo(items)
            self.makeNewWindow(aligned)
            self.alignments.append(list(items.values()))
        else:
            self.mainStatus.showMessage("Alignment already opened", msecs=5000)

    def makeNewWindow(self, ali):
        sub = views.MDISubWindow()
        sub.setAttribute(Qt.WA_DeleteOnClose, False)
        sub.setWidget(views.AlignSubWindow(ali))
        print(sub.widget())

        # Save the window in the dictionary
        print("Checking if title present")
        query = self.windows.get(str(sub.windowTitle()))
        print(query)
        print("Checking if window present")
        query2 = self.windows.get(sub)
        print(query2)
        if not query:
            print("adding new window to dictionary")
            self.windows[str(sub.windowTitle())] = sub.widget()
            print(self.windows)
        else:
            count = self.queryWindows(sub.windowTitle())
            sub.setWindowTitle(sub.windowTitle() + " " + count)
            self.windows[sub.windowTitle()] = sub.widget()

        # Show window in the view panel
        self.mdiArea.addSubWindow(sub)
        node = models.SeqNode(sub.windowTitle(), window=sub)
        self.projectRoot.appendRow(node)
        sub.show()
        print("~~~~~~~~~~~DONE~")

    def queryWindows(self, query):
        """ Quick counter to prevent duplicate naming of windows"""
        count = []
        print("Current keys: ", self.windows.keys())
        for key in self.windows.keys():
            print(key)
            if query in key:
                count.append(key)
        return str(len(count))

    def reopenWindow(self):
        """
        Checks to see if a window is open already.
        If it is not, reopens the window. If it is, gives focus.
        Also refreshes the title. May want to set up a separate listener for that.
        """
        item = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
        try:
            sub = item.text()
            sub.setWindowTitle(item.text())
            if not sub.isVisible():
                sub.show()
            self.mdiArea.setActiveSubWindow(sub)
        except:
            self.mainLogger.debug("Lost the window!")
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
        seq_GPI1B = Bseq.MutableSeq(test2[1],generic_protein)
        test2alt = [test2[0], seq_GPI1B]
        test = [test1alt, test2alt]
        for i in list(range(0, len(test))):
            node = models.SeqNode(test[i][0], sequence=test[i][1])
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
