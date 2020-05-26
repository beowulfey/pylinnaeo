import logging
import sys
import traceback

from PyQt5.QtCore import pyqtSignal, QThread, QTimer
import clustalo

"""
Additional classes and functions that are used within Sherlock, but are not responsible for viewing data.
"""

class SeqThread(QThread):
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
    """ Tool for checking a title list. Used for generating new titles if duplicate"""
    #print("\nBEGIN CHECK-- Layer ", layer)
    #print("Searching in", titles)
    if name not in titles:
        # SAFE! You can add and return
    #    print("Final:",name)
        finalname = name
    elif name[-2] == "_" and int(name[-1]):
        # if there's already a name with an _1, add a number
        newlayer = layer+1
        newname = str(name[:-1] + str(newlayer))
    #    print("Trying: ",newname)
        finalname, titles = checkName(newname, titles, layer=newlayer)
        if layer > 0:
    #        print("returning")
            return finalname, titles
    else:
        # It's a duplicate! Better
    #    print("Dupe found: [", name, "] Descending")
        newlayer = layer + 1
        newname = str(name + "_" + str(newlayer))
    #    print("Trying ", newname)
        # Run the check again with the new name
        finalname, titles = checkName(newname, titles, layer=newlayer)
        if layer > 0:
            return finalname, titles
    if layer == 0:
        titles.append(finalname)
    #    print("Appended: ",titles, "\n")
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


def buildRuler(chars, gap, start, end):
    """
    Fairly complicated function for generating the ruler. Calculates the spacing
    and labeling based on the width of the screen. Might be computational intensive,
    so I make a point to hide it (and the colors) when resizing.
    """
    ruler = None
    if start != end:
        #print("\nBEGIN LINE!")
        #print("Chars, Gap, Start, End: ", chars,gap,start,end)
        if gap != 0 and chars != end-start:
            if end - start > 8:
                chars = (end - start)
        #        print("Setting chars to ",chars)
            else:
                chars = 8
        labels = list(range(start+1, end+1))
        spacing = chars
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
        #print("Selected spacing: ", spacing)
        n = 0
        speclabels = []
        rulerlist = []
        while n < chars:
            speclabels.append(labels[n])
            n+=spacing
        if labels[-1] not in speclabels:
            speclabels.append(labels[-1])
        for n in range(1,4):
            if labels[-1]-n in speclabels:
                speclabels.remove(labels[-1]-n)
        # print(speclabels)
        if len(speclabels)>2:
            count = range(1,len(speclabels))
            # print("COUNT", list(count))
            for x in count:
                if x == count[-1]:
                    # print("Appending ", str(speclabels[x - 1]))
                    rulerlist.append(str(speclabels[x - 1]))
                    # print("FINAL gap, ", (chars - len(str("".join(rulerlist)))-len(str(speclabels[x]))))
                    rulerlist.append(" "*(chars - len(str("".join(rulerlist)))-len(str(speclabels[x]))))
                    # print("Appending final:", str(speclabels[x]))
                    rulerlist.append(str(speclabels[x]))
                else:
                    # print("Appending ", str(speclabels[x - 1]))
                    rulerlist.append(str(speclabels[x-1]))
                    first = 0
                    if x == count[0]:
                        first = len(str(speclabels[x-1]))
                    next = len(str(speclabels[x]))
                    # print("For gap after", str(speclabels[x - 1]), "space is",(spacing - next - first+1))
                    rulerlist.append(" "*(spacing - next - first+1))
            ruler = "".join(rulerlist)
        elif len(speclabels)==2:
            ruler = str(start + 1) + " " * (chars - len(str(start + 1)) - len(str(end))) + str(end)
        else:
            ruler = str(start+1)
    return ruler


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
            result = clustalo.clustalo(*self.args, **self.kwargs)
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
