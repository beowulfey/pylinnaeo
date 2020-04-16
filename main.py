#!/usr/bin/python3

from PyQt5.QtWidgets import QApplication, QTabWidget
import sys

import bioglot
import bioglotui

class BioGlot(QTabWidget, bioglotui.Ui_TranslateTest):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.translateButton.clicked.connect(self.translate_sequence)

    def translate_sequence(self):
        user_sequence = self.dnaSequence.toPlainText()
        contained_seq = bioglot.DnaContainer()
        contained_seq.setSequence(user_sequence)
        translated = contained_seq.getTranslation()
        self.proteinSequence.setPlainText(translated)




def main():
    app = QApplication(sys.argv)
    window = BioGlot()
    window.show()
    sys.exit(app.exec_())






if __name__ == '__main__':
    main()





