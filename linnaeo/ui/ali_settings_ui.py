# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui/ali_settings.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(213, 527)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(Form)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_3 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.line_2 = QtWidgets.QFrame(Form)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.verticalLayout.addWidget(self.line_2)
        self.checkRuler = QtWidgets.QCheckBox(Form)
        self.checkRuler.setObjectName("checkRuler")
        self.verticalLayout.addWidget(self.checkRuler)
        self.checkColors = QtWidgets.QCheckBox(Form)
        self.checkColors.setObjectName("checkColors")
        self.verticalLayout.addWidget(self.checkColors)
        self.label_2 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.line = QtWidgets.QFrame(Form)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout.addWidget(self.line)
        self.comboTheme = QtWidgets.QComboBox(Form)
        self.comboTheme.setObjectName("comboTheme")
        self.comboTheme.addItem("")
        self.comboTheme.addItem("")
        self.comboTheme.addItem("")
        self.comboTheme.addItem("")
        self.comboTheme.addItem("")
        self.comboTheme.addItem("")
        self.verticalLayout.addWidget(self.comboTheme)
        self.checkConsv = QtWidgets.QCheckBox(Form)
        self.checkConsv.setEnabled(False)
        self.checkConsv.setObjectName("checkConsv")
        self.verticalLayout.addWidget(self.checkConsv)
        self.checkColorDesc = QtWidgets.QCheckBox(Form)
        self.checkColorDesc.setEnabled(False)
        self.checkColorDesc.setObjectName("checkColorDesc")
        self.verticalLayout.addWidget(self.checkColorDesc)
        self.label = QtWidgets.QLabel(Form)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.comboFont = QtWidgets.QFontComboBox(Form)
        font = QtGui.QFont()
        font.setFamily("Liberation Mono")
        self.comboFont.setFont(font)
        self.comboFont.setWritingSystem(QtGui.QFontDatabase.Latin)
        self.comboFont.setFontFilters(QtWidgets.QFontComboBox.MonospacedFonts)
        font = QtGui.QFont()
        font.setFamily("Andale Mono")
        self.comboFont.setCurrentFont(font)
        self.comboFont.setObjectName("comboFont")
        self.verticalLayout.addWidget(self.comboFont)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.spinFontSize = QtWidgets.QSpinBox(Form)
        self.spinFontSize.setObjectName("spinFontSize")
        self.horizontalLayout.addWidget(self.spinFontSize)
        self.label_4 = QtWidgets.QLabel(Form)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout.addWidget(self.label_4)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout)
        spacerItem1 = QtWidgets.QSpacerItem(20, 200, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label_3.setText(_translate("Form", "Window Settings"))
        self.checkRuler.setText(_translate("Form", "Show ruler"))
        self.checkColors.setText(_translate("Form", "Show colors"))
        self.label_2.setText(_translate("Form", "Theme"))
        self.comboTheme.setCurrentText(_translate("Form", "Default"))
        self.comboTheme.setItemText(0, _translate("Form", "Default"))
        self.comboTheme.setItemText(1, _translate("Form", "Bold"))
        self.comboTheme.setItemText(2, _translate("Form", "ColorSafe"))
        self.comboTheme.setItemText(3, _translate("Form", "Monochrome"))
        self.comboTheme.setItemText(4, _translate("Form", "Grayscale"))
        self.comboTheme.setItemText(5, _translate("Form", "Rainbow"))
        self.checkConsv.setText(_translate("Form", "By conservation"))
        self.checkColorDesc.setText(_translate("Form", "Color meanings"))
        self.label.setText(_translate("Form", "Font:"))
        self.label_4.setText(_translate("Form", "pts"))
