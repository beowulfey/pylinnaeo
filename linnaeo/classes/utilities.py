import json
import logging
import re
import sys
import traceback
import urllib
from urllib.error import HTTPError

import clustalo
import numpy
from Bio import PDB
from Bio.PDB import DSSP, PDBIO
from PyQt5.QtCore import pyqtSignal, QThread, QTimer, Qt, QTemporaryFile
from bioservices import UniProt

"""
Additional classes and functions that are used within Linnaeo, but are not responsible for viewing data.
"""

def checkConservation(res, ref):
    """ Reads in a residue and checks whether it falls within a conserved category relative to the reference. """
    conserved = {
        # Adapted from http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
        # Original paper: https://pubmed.ncbi.nlm.nih.gov/14872534/
        # Direct conserved first --> exact category
        0: ['A','I','L','V','M'], 1: ['R', 'H', 'K'], 2: ['S', 'T'],
        3: ['D', 'E'], 4: ['N','Q'], 5: ['C', 'C'], 6: ['G','G'], 7: ['P','P'], 8: ['W','W'],
        9: ['Y','Y'], 10: ['F','F'],
        # Very similar; mostly compatible swaps
        11: ['I', 'L', 'V', 'M', 'F',], 12: ['A', 'S'], 13: ['D', 'N'], 14: ['E', 'Q'],
        15: ['F', 'W', 'Y'], 16: ['A', 'T'], 17: ['H', 'Q'],
        # Less similar; sort of compatible swaps
        18: ['R', 'E', 'Q'], 19: ['K','E','Q'], 20: ['C','S',], 21: ['T','I','L','V','M','F'],
        22: ['H', 'Y'], 23: ['W','Q']
    }
    for key, value in conserved.items():
        #print("Searching %s, %s, in %s" % (res, ref, value))
        if res in value and ref in value:
        #    print('%s matches %s within %s' % (res, ref, value))
            del value, res, ref, conserved
            return key
        else:
            continue
    del key, value, conserved, res, ref
    return None


def redrawBasic(seqs, chars, lines, rulers=False, dssp=False):
    """
    Black and white only; keeps space for ruler but does not calculate the ruler or DSSP stuff,
    which helps with speed during window size changes.
    """
    html = ["<pre>\n"]
    n = 0
    for line in range(lines):
        if rulers:
            html.append("\n")
        if dssp:
            html.append("\n")
        start = n * chars
        end = n * chars + chars
        if line == lines - 1:
            end = start + len(seqs[0][start:])
        for i in range(len(seqs)):
            html.append("".join([x[0][-8] for x in seqs[i][start:end]]))
            if line == lines-1:
                label = str(seqs[i][end - 1][1])
                if label == "0":
                    for y in range(chars):
                        label = str(seqs[i][end - 1 - y][1])
                        if label != "0":
                            break
                    del y
                html.append(" "*2 + label + "\n")
                del label
            else:
                html.append("\n")

        html.append("\n")
        n+=1
    html.append("</pre>")
    final = "".join(html)
    del html, seqs, chars, lines, rulers, dssp, i, line, start, end, n
    return final


def redrawFancy(seqs, chars, lines, rulers, colors, dssp):
    """
    Fancy like the name implies. Called at the end of resize events. Keeps your opinion on colors and rulers.
    Keeps a lot of whitespace because the tooltip and calculation of the residue ID depends on certain whitespace.
    """
    # All the possible symbols from DSSP. The unicode values were drawn by me in my modified default font.
    '''lookup = {'H': '&#x27B0;', 'G': '&nbsp;', 'I': '&nbsp;',
              '>': '&#x27B2;', 'E': '&#x27B1;', 'B': '&nbsp;',
              'T': '&#x27B3;', 'S': '&#x27B3;',
              '-': '&nbsp;', 'C': '&nbsp;'}'''
    # Only shows SS for classic helices and beta-strands for cleanliness
    lookup = {'H': '&#x27B0;', 'G': '&nbsp;', 'I': '&nbsp;',
              '>': '&#x27B2;', 'E': '&#x27B1;', 'B': '&nbsp;',
              'T': '&nbsp;', 'S': '&nbsp;',
              '-': '&nbsp;', 'C': '&nbsp;'}
    html = ["<pre>"]
    n = 0
    #print("Redraw Fancy: DSSP",dssp)
    for line in range(lines):
        start = n * chars
        end = start + len(seqs[0][start:]) if line == lines - 1 else n * chars + chars
        if rulers:
            # If there is a horizontal ruler, build the ruler and thread it in.
            html.append(str(buildRuler(chars, start, end))+"\n")
        for i in range(len(seqs)):
            if dssp and i == 0:
                ss = ['<span style=\"font-family:Default-Noto;font-size:inherit;\">']
                # TODO: THIS DOES NOT WORK FOR ALIGNMENTS; NEED TO LOOK AT THE TRUE INDEX NOT RAW INDEX NUMBER?!
                for index, x in enumerate(seqs[i][start:end]):
                    #print(index+start, x[2])
                    if x[2]:
                        if x[2] == "E":
                            #print(x[2],index+start,"Next residue:", seqs[i][index+start+1][2], index+start+1)
                            if seqs[i][index+start+1][2] != 'E':
                                ss.append(lookup['>'])
                            #    print('using arrowhead')
                            else:
                            #    print("using rectangle")
                                ss.append(lookup[x[2]])
                        else:
                            ss.append(lookup[x[2]])
                ss.append("</span>")
                html.append("".join(ss))
                #html.append("".join([x[2] if x[2] else "&nbsp;" for x in seqs[i][start:end]]))
                if line == lines-1:
                    html.append("&nbsp;"*(chars-(end-start)))
                html.append("\n")
            if colors:
                # The whole thing has the HTML color as well
                html.append("".join([x[0] for x in seqs[i][start:end]]))
            else:
                # Only use the residue, not the color HTML part.
                html.append("".join([x[0][-8] for x in seqs[i][start:end]]))
            if line == lines-1:
                # If last line, append the final residue ID number.
                label = str(seqs[i][end - 1][1])
                if label == "0":
                    for y in range(chars):
                        label = str(seqs[i][end - 1 - y][1])
                        if label != "0":
                            break
                html.append("&nbsp;"*2 + label + "&nbsp;"*(chars-(end-start)-(len(label)+2))+"\n")
            else:
                html.append("\n")
        if line < lines - 1:
            html.append(" "*chars+"\n")
        n += 1
    html.append("</pre>")
    final = "".join(html)
    del seqs, chars, lines, rulers, colors, dssp, html, line, start, end, i,
    return final


def buildRuler(chars, start, end):
    """ Semi-complicated function for building the horizontal ruler. Divides the screen width into even intervals. """
    ruler = None
    if start != end:
        interval = chars
        if 25 <= chars <= 50:
            interval = int(chars / 2)
        elif 50 < chars <= 100:
            interval = int(chars / 3)
        elif 100 < chars <= 150:
            interval = int(chars / 4)
        elif 150 < chars <= 250:
            interval = int(chars / 5)
        elif 250 < chars:
            interval = int(chars/6)
        labels = list(numpy.arange(start, end, interval))
        if len(labels) == 1:
            #    print(labels[0])
            if labels[0] <= (start + 1):
                #        print("removed ", labels[0])
                labels.pop(0)
            elif labels[0] - (start + 1) < len(str(labels[-1])) * 2 + 2:
                #        print("removed ", labels[0])
                labels.pop(0)
        if len(labels) > 1:
            rm = []
            for x in range(2):
                #        print(labels[x])
                if labels[x] - (start + 1) < len(str(labels[-1])) * 2 + 2:
                    #            print("Removed ",labels[x])
                    rm.append(labels[x])
                    if len(labels) <= 1:
                        break
                if end - labels[-x] < len(str(labels[-1])) * 2 + 2:
                    #            print("Removed ",labels[-x])
                    rm.append(labels[-x])
            [labels.remove(x) for x in rm]
            del rm
        labels.insert(0, start + 1)
        if end - (start + 1) >= len(str(start + 1)) + len(str(end)) + 2:
            labels.append(end)
        # print(labels)
        rulel = ['<u>' + str(labels[0])[0] + '</u>' + str(labels[0])[1:]]
        for x in labels[1:]:
            prevx = labels[labels.index(x) - 1]
            xlab = len(str(x))
            if labels.index(x) == 1 and len(str(prevx)) > 1:
                # Different for first space because I left align the first number.
                xlab = len(str(x)) + len(str(prevx)) - 1
            spacer = x - prevx - xlab
            #   print("Spacer is (%s-%s)-%s=%s " % (x, prevx, xlab, spacer))
            rulel.append("&nbsp;" * spacer)
            rulel.append(str(x)[:-1] + '<u>' + str(x)[-1:] + '</u>')
            del spacer
            del x, prevx, xlab,
        count = list(re.findall(r'(&nbsp;)|[0-9]',"".join(rulel)))
        rulel.append("&nbsp;"*(chars-len(count)))
        ruler = "".join(rulel)
        del rulel, count, labels, interval
    del chars, start, end,
    return ruler


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
    del layer, name
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


def nodeSelector(tree, model):
    """ Utility for selection of nodes and anything under a folder. """
    seqs = []
    copied = []
    indices = tree.selectedIndexes()
    if indices:
        for index in indices:
            node = model.itemFromIndex(index)
            if not node.data(role=Qt.UserRole + 2):
                for subnode in iterTreeView(model.itemFromIndex(index)):
                    i = model.indexFromItem(subnode)
                    copied.append(i)
                    seqr = subnode.data(role=Qt.UserRole + 2)[0]
                    seqs.append(seqr)
                    del seqr
            if node.data(role=Qt.UserRole + 2) and index not in copied:
                copied.append(index)
                seqr = node.data(role=Qt.UserRole + 2)[0]
                seqs.append(seqr)
                del seqr
            del node
        del index
    del tree, model, copied
    return indices, seqs


class SeqThread(QThread):
    """
     Determines which type of redraw should occur based on the values of rulers, colors, and "fancy", then
     does this in a separate thread. See redraw functions above.
     """
    finished = pyqtSignal()

    def __init__(self, *args, fancy=True, parent=None):
        QThread.__init__(self, parent)
        self.seqLogger = logging.getLogger("SEQDRAW")
        self.html = None
        self.seqs = args[0]
        self.chars = args[1]
        self.lines = args[2]
        self.rulers = args[3]
        self.colors = args[4]
        self.dssp = args[5]
        self.fancy = fancy
        #self.seqLogger.debug("Sequence thread created")

    def run(self):
        while not self.html:
            if self.fancy:
                self.html = redrawFancy(self.seqs, self.chars, self.lines, self.rulers, self.colors, self.dssp)
            else:
                self.html = redrawBasic(self.seqs, self.chars, self.lines, self.rulers, self.dssp)
        #self.seqLogger.debug("Sequence thread finished")
        self.finished.emit()



class GetPDBThread(QThread):
    """
    Input is a tuple of [SeqRs],Indices. To keep it compatible, if this is being called from a different
    method other than the button (like from , it will ignore the fact that
    """
    finished = pyqtSignal(list)
    logger = pyqtSignal(str)
    
    def __init__(self, input, parent=None):
        QThread.__init__(self, parent)
        self.seqs = input[0]
        self.nodes = None if not input[1] else input[1]
        self.u = UniProt(verbose=False)
        self.PDBLogger = logging.getLogger("PDBSearch")

    def run(self):
        returned = []
        # TODO: Should pop up modal if not confident in the structure results!
        index = None
        for seq in self.seqs:
            if self.nodes:
                index = self.nodes[self.seqs.index(seq)]
            pid = seq.id
            struct = None
            structseq = None
            self.PDBLogger.info("Searching with ID %s" % pid)
            self.logger.emit("Searching with ID %s" % pid)
            uniID = self.u.search(pid, columns="id, genes, organism, protein names")
            alignoff = 0
            if uniID:
                self.PDBLogger.info('Results!\n\n%s' % uniID)
                self.logger.emit('Results collected for search: %s' % uniID)
                result = uniID.split("\n")
                ids = []
                for line in result[1:]:
                    ids.append(line.split("\t"))
                coordinates = None
                i = 0
                while coordinates == None:
                    self.PDBLogger.info('Attempting search with %s from %s' %  (ids[i][1], ids[i][2]) )
                    self.logger.emit('Attempting search with %s from %s' %  (ids[i][1], ids[i][2]) )
                    structurl = "https://swissmodel.expasy.org/repository/uniprot/%s.json" % ids[i][0]
                    self.PDBLogger.debug('Searching SwissModel repository: %s' % structurl)
                    self.logger.emit('Searching SwissModel for structure')
                    try:
                        with urllib.request.urlopen(structurl) as url:
                            data = json.loads(url.read().decode())
                        if data['result']:
                            #print("Data found")
                            result = data['result']
                            if result['structures']:
                                #print("structures found")
                                structs = result['structures']
                                structseq = result['sequence']
                                self.PDBLogger.info("QUERY: \n%s" % str(seq.seq).replace("-",""))
                                self.PDBLogger.info("RESULT: \n%s" % structseq)
                                if str(seq.seq).replace("-", "") == structseq:
                                    # They should match, else keep looking
                                    if structs[0]:
                                        # print("accessing first model")
                                        struct0 = structs[0]
                                        if struct0['coordinates']:
                                            coordinates = struct0['coordinates']
                                            alignoff = int(struct0['chains']['A'][0]['uniprot']['from']) - 1
                                            self.PDBLogger.debug("MODEL ACQUIRED")
                                            self.logger.emit('MODEL ACQUIRED!')
                                        else:
                                            i += 1
                                            continue
                                    else:
                                        i += 1
                                        continue
                                else:
                                    self.PDBLogger.debug("Seq didn't match, trying with next model")
                                    i += 1
                                    continue
                            else:
                                i += 1
                                continue
                        elif i == len(ids):
                            self.PDBLogger.info("Sorry, no models found")
                            break
                        else:
                            i += 1
                            continue
                    except HTTPError:
                        break

                if coordinates:
                    offset = 0
                    start = structseq[:7]
                    #print(start)
                    for x in range(len(seq)):
                        end = x + 7
                        if str(seq.seq)[x:end] == start:
                            self.PDBLogger.debug("Sequence offset is %s residues" % x)
                            offset = x + alignoff
                            self.PDBLogger.info("Alignment offset is %s residues" % offset)
                            self.logger.emit("Alignment offset is %s residues" % offset)
                    parser = PDB.PDBParser()
                    tmp = QTemporaryFile()
                    with urllib.request.urlopen(coordinates) as url:
                        myfile = url.read()
                        if tmp.open():
                            tmp.write(myfile)
                            struct = parser.get_structure(ids[1], tmp.fileName())
                            self.PDBLogger.debug("STRUCTURE PARSED")
                            self.logger.emit("STRUCTURE PARSED")
                            #print(struct, type(struct))
                            returned.append([struct, seq, index, offset])

                else:
                    self.PDBLogger.debug("Sorry, no models found!!!")
            else:
                self.PDBLogger.info("NO STRUCTURE FOUND")
                self.logger.emit("No structure found, sorry!")
            self.finished.emit(returned)


class DSSPThread(QThread):
    """ Uses a stored PDB to calculate the secondary structure of a sequence using DSSP """
    finished = pyqtSignal(list)

    def __init__(self, *args, parent=None):
        QThread.__init__(self, parent)
        self.args = args
        self.struct = self.args[0][0]
        self.seq = self.args[0][1]
        self.node = self.args[0][2]
        self.offset = self.args[0][3]
        self.result = {}

    def run(self):
        tmp = QTemporaryFile()
        result = {}
        if tmp.open():
            io = PDBIO()
            io.set_structure(self.struct)
            io.save(tmp.fileName())
            dssp = DSSP(self.struct[0], tmp.fileName(), dssp='mkdssp')
            prevChain = next(iter(dssp.keys()))[0]
            for key in dssp.keys():
                #print(key[0])
                if key[0] ==prevChain:
                    #print(key)
                # I THINK I'M DOING THIS PART WRONG
                    result[dssp[key][0]+self.offset] = dssp[key][2]
        self.finished.emit([result, self.seq, self.node])


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
        self.args = list(args)
        parent = self.args[0]
        self.args.pop(0)
        QThread.__init__(self, parent)
        self.kwargs = kwargs
        self.aligned = {}
        del parent
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
            self.clustalLogger.debug("Thread returned alignment successfully")


class ProcTimerThread(QThread):
    """
    Thread for the timer, because Windows complains like hell otherwise.
    For memory and CPU usage refresh.
    """
    timeout = pyqtSignal()

    def __init__(self, parent=None):
        QThread.__init__(self, parent)
        self.processTimer = QTimer()
        self.processTimer.setInterval(1000)
        self.processTimer.timeout.connect(self.timerDone)
        self.processTimer.start()

    def timerDone(self):
        self.timeout.emit()



