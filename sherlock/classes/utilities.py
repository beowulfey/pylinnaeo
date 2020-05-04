import logging
import sys
import traceback
from PyQt5.QtCore import QObject, pyqtSignal, QRunnable, pyqtSlot, QThread
from clustalo import clustalo

"""
Additional classes and functions that are used within Sherlock, but are not responsible for viewing data.
"""


def checkName(name, titles, layer=0):
    """ Tool for checking a title list. Used for generating new titles if duplicate"""
    print(name)
    print(titles)
    if name not in titles:
        print("SAFE!")
        # SAFE! You can add and return
        finalname = name
        titles.append(finalname)
    elif name[-2] == "_" and int(name[-1]):
        print("ADDING NUMBER")
        # if there's already a name with an _1, add a number
        finalname = str(name[:-1] + str(int(name[-1]) + 1))
        titles.append(finalname)
    else:
        print("IS DUPLICATE! Check new name")
        # It's a duplicate! Better
        finalname = str(name + "_" + str(titles.count(name)))
        newlayer = layer + 1
        # Run the check again with the new name
        finalname, titles = checkName(finalname, titles, layer=newlayer)
        titles.append(finalname)
    return finalname, titles


def iterTreeView(root):
    """
    Internal function for iterating a TreeModel.
    Usage: for node in _iterTreeView(root): etc.
    """

    def recurse(parent):
        for row in range(parent.rowCount()):
            child = parent.child(row)
            yield child
            if child.hasChildren():
                yield from recurse(child)

    if root is not None:
        yield from recurse(root)


class AlignThread(QThread):
    finished = pyqtSignal()
    error = pyqtSignal(tuple)

    def __init__(self, *args, **kwargs):
        self.clustalLogger = logging.getLogger("ClustalO")
        QThread.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self.aligned = {}
        #self.clustalLogger.debug("Thread for ClustalO created")

    def run(self):
        try:
            result = clustalo(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.aligned = result
            self.clustalLogger.debug("ClustalO thread returned alignment successfully")


