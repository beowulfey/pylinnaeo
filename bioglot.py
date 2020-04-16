#!/usr/bin/python3

# BIOGLOT
# This is my shitty code for managing all the parts I'm going to somehow
# hack together into a piece of software, maybe.

from PyQt5.QtWidgets import QMainWindow
from ui import bioglot_ui
import biotite.sequence


class BioGlot(QMainWindow, bioglot_ui.Ui_MainWindow):
    # Intitalize UI
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)

    # UI interaction code goes here
    

#    def translate_sequence(self):
#        user_sequence = self.dnaSequence.toPlainText()
#        contained_seq = bioglot.DnaContainer()
#        contained_seq.setSequence(user_sequence)
#        translated = contained_seq.getTranslation()
#        self.proteinSequence.setPlainText(translated)
