# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'quit.ui'
#
# Created by: PyQt5 UI code generator 5.14.2
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_closeConfirm(object):
    def setupUi(self, closeConfirm):
        closeConfirm.setObjectName("closeConfirm")
        closeConfirm.setWindowModality(QtCore.Qt.ApplicationModal)
        closeConfirm.resize(300, 100)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(closeConfirm.sizePolicy().hasHeightForWidth())
        closeConfirm.setSizePolicy(sizePolicy)
        self.gridLayout = QtWidgets.QGridLayout(closeConfirm)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(closeConfirm)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(closeConfirm)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Discard|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 1, 0, 1, 1)

        self.retranslateUi(closeConfirm)
        self.buttonBox.accepted.connect(closeConfirm.accept)
        self.buttonBox.rejected.connect(closeConfirm.reject)
        QtCore.QMetaObject.connectSlotsByName(closeConfirm)

    def retranslateUi(self, closeConfirm):
        _translate = QtCore.QCoreApplication.translate
        closeConfirm.setWindowTitle(_translate("closeConfirm", "Dialog"))
        self.label.setText(_translate("closeConfirm", "Save your workspace before closing?"))
