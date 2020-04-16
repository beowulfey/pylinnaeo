# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'bioglot.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_TranslateTest(object):
    def setupUi(self, TranslateTest):
        TranslateTest.setObjectName("TranslateTest")
        TranslateTest.resize(400, 300)
        self.dna = QtWidgets.QWidget()
        self.dna.setObjectName("dna")
        self.splitter = QtWidgets.QSplitter(self.dna)
        self.splitter.setGeometry(QtCore.QRect(70, 10, 256, 224))
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.splitter.setObjectName("splitter")
        self.dnaSequence = QtWidgets.QPlainTextEdit(self.splitter)
        self.dnaSequence.setTabChangesFocus(True)
        self.dnaSequence.setObjectName("dnaSequence")
        self.translateButton = QtWidgets.QPushButton(self.splitter)
        self.translateButton.setObjectName("translateButton")
        TranslateTest.addTab(self.dna, "")
        self.protein = QtWidgets.QWidget()
        self.protein.setObjectName("protein")
        self.proteinSequence = QtWidgets.QPlainTextEdit(self.protein)
        self.proteinSequence.setGeometry(QtCore.QRect(70, 10, 256, 192))
        self.proteinSequence.setReadOnly(True)
        self.proteinSequence.setCursorWidth(0)
        self.proteinSequence.setObjectName("proteinSequence")
        TranslateTest.addTab(self.protein, "")

        self.retranslateUi(TranslateTest)
        TranslateTest.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(TranslateTest)

    def retranslateUi(self, TranslateTest):
        _translate = QtCore.QCoreApplication.translate
        TranslateTest.setWindowTitle(_translate("TranslateTest", "TranslateTest"))
        self.dnaSequence.setPlaceholderText(_translate("TranslateTest", "Paste DNA sequence here"))
        self.translateButton.setText(_translate("TranslateTest", "Translate!"))
        TranslateTest.setTabText(TranslateTest.indexOf(self.dna), _translate("TranslateTest", "DNA"))
        TranslateTest.setTabText(TranslateTest.indexOf(self.protein), _translate("TranslateTest", "Protein"))

