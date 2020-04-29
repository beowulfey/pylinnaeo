from PyQt5 import QtCore, QtGui, QtWidgets

############################################
''' Classes '''
############################################
class Main_Window( QtWidgets.QDialog ):
    def __init__( self ):
        super( Main_Window, self ).__init__( )

        self.create_gui()
        self.create_layout()
        self.create_connections()
        self.get_contents()

    #--------------------------------------------------------------------
    def create_gui( self ):
        self.tv_model=MyModel()
        self.tv_file_list = File_List( self )

    #--------------------------------------------------------------------
    def create_layout( self ):
        self.main_layout = QtGui.QVBoxLayout( self )
        self.main_layout.addWidget( self.tv_file_list )
        self.setLayout( self.main_layout )

    #--------------------------------------------------------------------
    def get_contents(self):
        self.tv_model.clear()
        self.tv_model.setHorizontalHeaderLabels(["name","date"])
        contents=["path1","path2"]
        for path in contents:
            date = self.get_date(path)
            self.add_file(path,date)

    #--------------------------------------------------------------------
    def add_file(self, name, date):
        item1 = QtGui.QStandardItem(name)
        item2 = QtGui.QStandardItem(date)
        self.tv_model.appendRow([item1, item2])

    #--------------------------------------------------------------------
    def get_date(self, path):
        return "a date"

    #--------------------------------------------------------------------
    def create_connections( self ):
        self.tv_file_list.clicked.connect( self.on_click )

    # slots --------------------------------------------------------------
    #def on_click(self, item ):
    #    # print item from first column
    #    index = self.tv_file_list.selectedIndexes()[0]
    #    item = self.tv_model.itemFromIndex(index)
    #    print item

    def on_click(self, item ):
        # print item from first column
        index = self.tv_file_list.selectedIndexes()[0]
        item = self.tv_model.itemFromIndex(index).text()
        print(item)



############################################
class MyModel(QtGui.QStandardItemModel):
    def __init__(self, parent=None):
        super(MyModel, self).__init__(parent)
    #--------------------------------------------------------------------
    def flags(self, index):
        flag = QtCore.Qt.ItemIsEnabled
        if index.isValid():
            flag |= QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsUserCheckable
        return flag

############################################
class File_List( QtWidgets.QTreeView ):
    ''' Create the file filters '''
    def __init__( self, mainUIWindow ):
        super( File_List, self ).__init__( parent )

        self.setModel(mainUIWindow.tv_model)
        self.setIndentation(0)

############################################
if __name__ == "__main__":
    # workaround for a bug in maya
    try:
        tree_view_ui.close()
        tree_view_ui.deleteLater()
    except:
        pass

    tree_view_ui = Main_Window()
    tree_view_ui.show()

    try:
        tree_view_ui.show()
    except:
        tree_view_ui.close()
        tree_view_ui.deleteLater()
