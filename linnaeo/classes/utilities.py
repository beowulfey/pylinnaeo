import logging
import sys
import traceback
from textwrap import TextWrapper

from PyQt5.QtCore import pyqtSignal, QThread, QTimer

from clustalo import clustalo

"""
Additional classes and functions that are used within Sherlock, but are not responsible for viewing data.
"""


class SeqWrap(TextWrapper):
    def __init__(self, **kwargs):
        super().__init__(kwargs)

    def _handle_long_word(self, reversed_chunks, cur_line, cur_len, width):
        ### SUBCLASSED OVERRIDDEN
        """_handle_long_word(chunks : [string],
                             cur_line : [string],
                             cur_len : int, width : int)

        Handle a chunk of text (most likely a word, not whitespace) that
        is too long to fit in any line.
        """
        # Figure out when indent is larger than the specified width, and make
        # sure at least one character is stripped off on every pass
        if width < 1:
            space_left = 1
        else:
            space_left = width - cur_len
            print("SPACE LEFT: ",space_left)

        # If we're allowed to break long words, then do so: put as much
        # of the next chunk onto the current line as will fit.
        if self.break_long_words:
            cur_line.append(reversed_chunks[-1][:space_left])
            print(cur_line)
            reversed_chunks[-1] = reversed_chunks[-1][space_left:]

        # Otherwise, we have to preserve the long word intact.  Only add
        # it to the current line if there's nothing already there --
        # that minimizes how much we violate the width constraint.
        elif not cur_line:
            cur_line.append(reversed_chunks.pop())

        # If we're not allowed to break long words, and there's already
        # text on the current line, do nothing.  Next time through the
        # main loop of _wrap_chunks(), we'll wind up here again, but
        # cur_len will be zero, so the next line will be entirely
        # devoted to the long word that we can't handle right now.


def checkName(name, titles, layer=0):
    """ Tool for checking a title list. Used for generating new titles if duplicate"""
    print("\nBEGIN CHECK-- Layer ", layer)
    print("Searching in", titles)
    if name not in titles:
        # SAFE! You can add and return
        print("Final:",name)
        finalname = name
    elif name[-2] == "_" and int(name[-1]):
        # if there's already a name with an _1, add a number
        newlayer = layer+1
        newname = str(name[:-1] + str(newlayer))
        print("Trying: ",newname)
        finalname, titles = checkName(newname, titles, layer=newlayer)
        if layer > 0:
            print("returning")
            return finalname, titles
    else:
        # It's a duplicate! Better
        print("Dupe found: [", name, "] Descending")
        newlayer = layer + 1
        newname = str(name + "_" + str(newlayer))
        print("Trying ", newname)
        # Run the check again with the new name
        finalname, titles = checkName(newname, titles, layer=newlayer)
        if layer > 0:
            return finalname, titles
    if layer == 0:
        titles.append(finalname)
        print("Appended: ",titles, "\n")
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


class TimerThread(QThread):
    timeout = pyqtSignal()

    def __init__(self):
        QThread.__init__(self)
        self.processTimer = QTimer()
        self.processTimer.setInterval(1000)
        self.processTimer.start()
        self.processTimer.timeout.connect(self.timerDone)

    def timerDone(self):
        self.timeout.emit()
