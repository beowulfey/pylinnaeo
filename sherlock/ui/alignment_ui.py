# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'alignment.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_aliWindow(object):
    def setupUi(self, aliWindow):
        aliWindow.setObjectName("aliWindow")
        aliWindow.resize(500, 300)
        aliWindow.setAcceptDrops(True)
        aliWindow.setStyleSheet("")
        self.gridLayout = QtWidgets.QGridLayout(aliWindow)
        self.gridLayout.setContentsMargins(0, 2, 0, 0)
        self.gridLayout.setVerticalSpacing(3)
        self.gridLayout.setObjectName("gridLayout")
        self.line = QtWidgets.QFrame(aliWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.line.sizePolicy().hasHeightForWidth())
        self.line.setSizePolicy(sizePolicy)
        self.line.setFrameShadow(QtWidgets.QFrame.Raised)
        self.line.setLineWidth(1)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setObjectName("line")
        self.gridLayout.addWidget(self.line, 0, 0, 1, 2)
        self.namePane = QtWidgets.QTextEdit(aliWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.namePane.sizePolicy().hasHeightForWidth())
        self.namePane.setSizePolicy(sizePolicy)
        self.namePane.setMinimumSize(QtCore.QSize(30, 100))
        self.namePane.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.namePane.setBaseSize(QtCore.QSize(75, 0))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(238, 238, 236))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(238, 238, 236))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(239, 239, 239))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        self.namePane.setPalette(palette)
        font = QtGui.QFont()
        font.setFamily("Liberation Mono")
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.namePane.setFont(font)
        self.namePane.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.namePane.setStyleSheet("font: 10pt \"Liberation Mono\";")
        self.namePane.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.namePane.setFrameShadow(QtWidgets.QFrame.Raised)
        self.namePane.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.namePane.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.namePane.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.namePane.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.namePane.setReadOnly(True)
        self.namePane.setObjectName("namePane")
        self.gridLayout.addWidget(self.namePane, 1, 0, 1, 1)
        self.alignPane = QtWidgets.QTextEdit(aliWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.alignPane.sizePolicy().hasHeightForWidth())
        self.alignPane.setSizePolicy(sizePolicy)
        self.alignPane.setMinimumSize(QtCore.QSize(200, 100))
        font = QtGui.QFont()
        font.setFamily("Liberation Mono")
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.alignPane.setFont(font)
        self.alignPane.setStyleSheet("font: 10pt \"Liberation Mono\";")
        self.alignPane.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.alignPane.setFrameShadow(QtWidgets.QFrame.Raised)
        self.alignPane.setLineWidth(2)
        self.alignPane.setMidLineWidth(2)
        self.alignPane.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.alignPane.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.alignPane.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.alignPane.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.alignPane.setReadOnly(True)
        self.alignPane.setTextInteractionFlags(QtCore.Qt.TextSelectableByKeyboard|QtCore.Qt.TextSelectableByMouse)
        self.alignPane.setObjectName("alignPane")
        self.gridLayout.addWidget(self.alignPane, 1, 1, 1, 1)

        self.retranslateUi(aliWindow)
        QtCore.QMetaObject.connectSlotsByName(aliWindow)

    def retranslateUi(self, aliWindow):
        _translate = QtCore.QCoreApplication.translate
        aliWindow.setWindowTitle(_translate("aliWindow", "New Alignment"))
        self.namePane.setHtml(_translate("aliWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Liberation Mono\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p dir=\'rtl\' style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
