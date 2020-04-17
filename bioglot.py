#!/usr/bin/python3

# BIOGLOT
# This is my shitty code for managing all the parts I'm going to somehow
# hack together into a piece of software, maybe.

from PyQt5.QtCore import QFile, QIODevice
from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QMainWindow
from ui import bioglot_ui
import biotite.sequence
import workspace
import configparser
from pathlib import Path


class BioGlot(QMainWindow, bioglot_ui.Ui_MainWindow):
    # Intitalize UI
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        config = configparser.ConfigParser()

        # Read in config (linux)
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
            test = ['GPI1B', 'MIFELFRFIFRKKKMLGYLSDLIGTLFIGDSTEKAMSLSQDATFVELKRHVEANE'+
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
            # NEED TO REMAKE STANDARDITEMMODEL so it has two options.
            # Messing around with building a Tree Model.
            model = QStandardItemModel()
            root = model.invisibleRootItem()
            for i in list(range(0, 4)):
                node = QStandardItem(test[0])
                if i == 2:
                    node2 = QStandardItem(test[0])
                    root.appendRow(node2)
                    parent = root.child(i-1)
                    parent.appendRow(node)
                else:
                    root.appendRow(node)
            self.workspaceTree.setModel(model)

    # UI interaction code goes here

#    def translate_sequence(self):
#        user_sequence = self.dnaSequence.toPlainText()
#        contained_seq = bioglot.DnaContainer()
#        contained_seq.setSequence(user_sequence)
#        translated = contained_seq.getTranslation()
#        self.proteinSequence.setPlainText(translated)
