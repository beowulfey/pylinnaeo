import logging
import sys
import traceback

from PyQt5.QtCore import pyqtSignal, QThread, QTimer
import clustalo

"""
Additional classes and functions that are used within Linnaeo, but are not responsible for viewing data.
"""


class SeqThread(QThread):
    """
     Determines which type of redraw should occur based on the values of rulers, colors, and "fancy", then
     does this in a separate thread.
     """
    finished = pyqtSignal()

    def __init__(self, *args, fancy=True):
        QThread.__init__(self)
        self.seqLogger = logging.getLogger("SEQDRAW")
        self.html = None
        self.seqs = args[0]
        self.chars = args[1]
        self.lines = args[2]
        self.rulers = args[3]
        self.colors = args[4]
        self.fancy = fancy

    def run(self):
        if self.fancy:
            self.html = redrawFancy(self.seqs, self.chars, self.lines, self.rulers, self.colors)
        else:
            self.html = redrawBasic(self.seqs, self.chars, self.lines, self.rulers)


def redrawBasic(seqs, chars, lines, rulers=False):
    """ Black and white only; keeps space for ruler but does not calculate the ruler, which helps with speed."""
    html = ["<pre>\n"]
    n = 0
    for line in range(lines):
        if rulers:
            html.append("\n")
        start = n * chars
        end = n * chars + chars
        if line == lines - 1:
            end = start + len(seqs[0][start:])
        for i in range(len(seqs)):
            html.append("".join([x[0][-8] for x in seqs[i][start:end]])+"\n")
        html.append("\n")
        n+=1
    html.append("</pre>")
    return "".join(html)


def redrawFancy(seqs, chars, lines, rulers=True, colors=True):
    """ Fancy like the name implies. Called at the end of resize events. Keeps your opinion on colors and rulers. """
    html = ["<pre>"]
    n = 0
    for line in range(lines):
        start = n * chars
        end = start + len(seqs[0][start:]) if line == lines - 1 else n * chars + chars
        gap = (n * chars + chars) - end
        if rulers:
            html.append(str(buildRuler(chars, gap, start, end))+"\n")
        for i in range(len(seqs)):
            if colors:
                html.append("".join([x[0] for x in seqs[i][start:end]])+"\n")
            else:
                html.append("".join([x[0][-8] for x in seqs[i][start:end]]) + "\n")
        html.append("\n")
        n += 1
    html.append("</pre>")
    return "".join(html)


def checkName(name, titles, layer=0):
    """ Tool for checking a list of titles. Used for generating a new title if it is duplicated"""
    # TODO: This isn't perfect, as it prevents duplicates of folders too... but that won't be fixed here.
    # TODO: Fix it in the item model.
    if name not in titles:
        # SAFE! You can add name and return
        finalname = name
    elif name[-2] == "_" and int(name[-1]):
        # if there's already a name with an _1, add a number
        newlayer = layer+1
        newname = str(name[:-1] + str(newlayer))
        finalname, titles = checkName(newname, titles, layer=newlayer)
        if layer > 0:
            return finalname, titles
    else:
        # It's a duplicate! Give it an underscore.
        newlayer = layer + 1
        newname = str(name + "_" + str(newlayer))
        # Run the check again with the new name
        finalname, titles = checkName(newname, titles, layer=newlayer)
        if layer > 0:
            return finalname, titles
    if layer == 0:
        titles.append(finalname)
    return finalname, titles


def iterTreeView(root):
    """
    Internal function for iterating a TreeModel. Shamelessly stolen from StackOverflow. It is wonderful though.
    Usage: for node in iterTreeView(root): --> returns all the nodes.
    """
    def recurse(parent):
        for row in range(parent.rowCount()):
            child = parent.child(row)
            yield child
            if child.hasChildren():
                yield from recurse(child)
    if root is not None:
        yield from recurse(root)


def buildRuler(chars, gap, start, end):
    """
    Fairly complicated function for generating the ruler. Calculates the spacing
    and labeling based on the width of the screen. Pretty computationally intensive,
    so I make a point to hide it (and the colors) when resizing.
    """
    ruler = None
    if start != end:
        if gap != 0 and chars != end-start:
            # Have to adjust the spacing for the last line
            # Don't go below 8 chars so the numbers don't merge.
            # May want to extend to 10... rare but some seqs are >1000
            if end - start > 8:
                chars = (end - start)
            else:
                chars = 8
        # labels is all the possible numbers between the start and end
        labels = list(range(start+1, end+1))
        spacing = chars
        # My hacky way to add more numbers as the screen increases in size.
        # Spacing is the distance between numbers of the rulers.
        # Not all numbers are divisible by 2, 3, etc so there are uneven spaces, which I account for badly below.
        # TODO: Consider changing this to floor?
        if chars < 20:
            pass
        elif 20 <= chars < 60:
            spacing = int(chars/2)
        elif 60 <= chars < 100:
            spacing = int(chars/3)
        elif 100 <= chars < 150:
            spacing = int(chars/4)
        elif 150 <= chars:
            spacing = int(chars/5)
        n = 0
        speclabels = []
        rulerlist = []
        # spec labels are the numbers that are actually used based on the spacing.
        while n < chars:
            speclabels.append(labels[n])
            n+=spacing
        if labels[-1] not in speclabels:
            # a little hack because sometimes the end doesn't  get added
            speclabels.append(labels[-1])
        for n in range(1,4):
            # another little hack for when math makes it so there are two labels right next to each other at the end
            if labels[-1]-n in speclabels:
                speclabels.remove(labels[-1]-n)
        if len(speclabels)>2:
            # more complicated logic for when there are multiple labels
            # builds a list (compressed to a string) of the labels and spacing, accounting for the length of the numbers
            count = range(1,len(speclabels))
            for x in count:
                if x == count[-1]:
                    rulerlist.append(str(speclabels[x - 1]))
                    rulerlist.append(" "*(chars - len(str("".join(rulerlist)))-len(str(speclabels[x]))))
                    rulerlist.append(str(speclabels[x]))
                else:
                    rulerlist.append(str(speclabels[x-1]))
                    first = 0
                    if x == count[0]:
                        first = len(str(speclabels[x-1]))
                    next = len(str(speclabels[x]))
                    rulerlist.append(" "*(spacing - next - first+1))
            ruler = "".join(rulerlist)
        elif len(speclabels)==2:
            # Ah! so easy.
            ruler = str(start + 1) + " " * (chars - len(str(start + 1)) - len(str(end))) + str(end)
        else:
            # Just show a single number.
            ruler = str(start+1)
    return ruler


class AlignThread(QThread):
    """
    Clustal Omega is run in a separate thread. Currently have no idea how to access the alignment order;
    I'm hoping there is a way, rather than returning it with the input order.
    I can't find anything in the source code of ClustalO for the API though, sadly.
    """
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
            result = clustalo.clustalo(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.aligned = result
            self.clustalLogger.debug("ClustalO thread returned alignment successfully")


class ProcTimerThread(QThread):
    """
    Thread for the timer, because Windows complains like hell otherwise.
    For memory and CPU usage refresh.
    """
    timeout = pyqtSignal()

    def __init__(self):
        QThread.__init__(self)
        self.processTimer = QTimer()
        self.processTimer.setInterval(1000)
        self.processTimer.timeout.connect(self.timerDone)
        self.processTimer.start()

    def timerDone(self):
        self.timeout.emit()

class ResizeTimerThread(QThread):
    """
    Thread for the timer, because Windows complains like hell otherwise.
    Also used for resizing with a super short countdown
    """
    timeout = pyqtSignal()

    def __init__(self):
        QThread.__init__(self)
        self.processTimer = QTimer()
        #self.processTimer.setSingleShot(True)
        self.processTimer.setInterval(500)
        self.processTimer.timeout.connect(self.timerDone)

    def timerDone(self):
        print(self.processTimer.interval(), "TIMER STOP")
        self.timeout.emit()

    def run(self):
        print("STARTING TIMER")
        self.processTimer.start()



