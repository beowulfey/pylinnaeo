import sys
import traceback
from PyQt5.QtCore import QObject, pyqtSignal, QRunnable, pyqtSlot, QThread
from clustalo import clustalo

"""
Additional classes and functions that are used within Sherlock, but are not responsible for viewing data.
"""


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
        QThread.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self.aligned = {}
        print("Thread created")

    def run(self):
        try:
            result = clustalo(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.aligned = result


