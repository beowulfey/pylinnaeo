import logging
import sys
import traceback
from textwrap import TextWrapper
import re

from PyQt5.QtCore import pyqtSignal, QThread, QTimer, QObject

from clustalo import clustalo

"""
Additional classes and functions that are used within Sherlock, but are not responsible for viewing data.
"""

"""
class SeqWrap(TextWrapper):
    _whitespace = '\t\n\x0b\x0c\r '
    word_punct = r'[\w!"\'&.,?]'
    letter = r'[^\d\W]'
    whitespace = r'[%s]' % re.escape(_whitespace)
    nowhitespace = '[^' + whitespace[1:]
    wordsep_re = re.compile(fr'''
        ( # any whitespace
          {whitespace}+
        #| # em-dash between words
        #  (?<={word_punct}) -{{2,}} (?=\w)
        | # word, possibly hyphenated
          {nowhitespace}+? #(?:
            # hyphenated word
            #  -(?: (?<={letter}{{2}}-) | (?<={letter}-{letter}-))
            #  (?= {letter} -? {letter})
            #| # end of word
            #  (?={whitespace}|\Z)
            #| # em-dash
            #  (?<={word_punct}) (?=-{{2,}}\w)
            #)
        )''',
                            re.VERBOSE)
    del word_punct, letter, nowhitespace

    def __init__(self, **kwargs):
        super().__init__(kwargs)

    def _handle_long_word(self, reversed_chunks, cur_line, cur_len, width):
        # OVERRIDDEN
        if width < 1:
            space_left = 1
        else:
            space_left = width - cur_len
        if self.break_long_words:
            cur_line.append(reversed_chunks[-1][:space_left])
            reversed_chunks[-1] = reversed_chunks[-1][space_left:]

        elif not cur_line:
            cur_line.append(reversed_chunks.pop())

"""


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
    ruler = None
    if start != end:
        print("\nBEGIN LINE!")
        print("Chars, Gap, Start, End: ", chars,gap,start,end)
        if gap != 0 and chars != end-start:
            if end - start > 8:
                chars = (end - start)
                print("Setting chars to ",chars)
            else:
                chars = 8
        labels = list(range(start+1, end+1))
        print(labels)
        #print(len(labels))
        # This shows numbers on beginning and end. Easy!
        # This shows them spaced out a set number.
        # Finds all the factors of the spacing.
        #factors = []
        #for i in range(1, chars+1):
        #    if chars % i == 0:
        #        factors.append(i)
        #print("Factors: ",factors)
        #pos = int(len(factors)/2)

        if chars < 20:
            spacing = chars
        elif 20 <= chars < 60:
            spacing = int(chars/2)
        elif 60 <= chars < 100:
            spacing = int(chars/3)
        elif 100 <= chars < 150:
            spacing = int(chars/4)
        elif 150 <= chars:
            spacing = int(chars/5)
        print("Selected spacing: ", spacing)
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
        print(speclabels)
        if len(speclabels)>2:
            count = range(1,len(speclabels))
            print("COUNT", list(count))
            for x in count:
                if x == count[-1]:
                    print("Appending ", str(speclabels[x - 1]))
                    rulerlist.append(str(speclabels[x - 1]))
                    print("FINAL gap, ", (chars - len(str("".join(rulerlist)))-len(str(speclabels[x]))))
                    rulerlist.append(" "*(chars - len(str("".join(rulerlist)))-len(str(speclabels[x]))))
                    print("Appending final:", str(speclabels[x]))
                    rulerlist.append(str(speclabels[x]))

                else:
                    print("Appending ", str(speclabels[x - 1]))
                    rulerlist.append(str(speclabels[x-1]))
                    first = 0
                    if x == count[0]:
                        first = len(str(speclabels[x-1]))
                    next = len(str(speclabels[x]))
                    print("For gap after", str(speclabels[x - 1]), "space is",(spacing - next - first+1))
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
