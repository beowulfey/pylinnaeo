import traceback
import copy
import time
from io import StringIO

import Bio
from Bio import SeqIO, AlignIO
from Bio.Alphabet import generic_protein
from Bio.Seq import MutableSeq, Seq
from Bio.Blast import NCBIWWW, NCBIXML
from PyQt5.QtCore import QFile, QIODevice, QDataStream, Qt, QDir, QTimer
from PyQt5.QtGui import QStandardItem
from PyQt5.QtWidgets import QFileDialog, QApplication
from bioservices import NCBIblast

from linnaeo.classes import widgets, models, utilities
from linnaeo.classes.displays import AboutDialog


class Slots:
    """
    Split out methods for clarity. These are all simple slot functions associated with UI actions!
    """

    def increaseTextSize(self):
        self._currentWindow.widget().increaseFont()
        self.mainStatus.showMessage("Setting font size to %s " % self._currentWindow.widget().getFontSize(), msecs=1000)

    def decreaseTextSize(self):
        self._currentWindow.widget().decreaseFont()
        self.mainStatus.showMessage("Setting font size to %s " % self._currentWindow.widget().getFontSize(), msecs=1000)

    def setCurrentWindow(self):
        self._currentWindow = self.mdiArea.activeSubWindow()
        if self._currentWindow:
            if self._currentWindow.widget().dssps:
                self.optionsPane.structureActivate(True)
            elif not self._currentWindow.widget().dssps:
                self.optionsPane.structureActivate(False)

    def showAbout(self):
        qDialog = AboutDialog(self)
        qDialog.exec()

    def saveImage(self):
        if self._currentWindow:
            self._currentWindow.widget().grab()  # TODO The first time this runs it redraws the window, but never after...
            w = self._currentWindow.widget().alignPane.size().width()
            sw = self._currentWindow.widget().alignPane.verticalScrollBar().size().width()
            self._currentWindow.widget().alignPane.verticalScrollBar().hide()
            self._currentWindow.widget().alignPane.size().setWidth(int(w - sw))
            file = QFileDialog.getSaveFileName(self, "Save as...", "name",
                                               "PNG (*.png);; BMP (*.bmp);;TIFF (*.tiff *.tif);; JPEG (*.jpg *.jpeg)")
            self._currentWindow.widget().grab().save(file[0] + file[1][-5:-1])
            self._currentWindow.widget().alignPane.verticalScrollBar().show()

    def toggleOptionsPane(self, state):
        self.optionsPane.show() if state else self.optionsPane.hide()

    def toggleRuler(self, state):
        if self._currentWindow:
            self._currentWindow.widget().toggleRuler(state)

    def toggleColors(self, state):
        if self._currentWindow:
            self._currentWindow.widget().toggleColors(state)

    def toggleStructure(self, state):
        print('Structure toggled, sending to window')
        if self._currentWindow:
            print("UPDATING STRUCTURE OF WINDOW", state)
            self._currentWindow.widget().toggleStructure(bool(state))

    def changeTheme(self):
        #try
        if self._currentWindow:
            self._currentWindow.widget().setTheme(self.optionsPane.comboTheme.currentText())
        #except:
        #    print("Skipping theme set")

    def changeFont(self, font):
        #try:
        if self._currentWindow:
            self._currentWindow.widget().setFont(font)
        #except:
        #    print("Skipping font set")

    def changeFontSize(self, size):
        if self._currentWindow:
            self._currentWindow.widget().setFontSize(size)

    def refreshParams(self, window):
        window.widget().setParams(self.optionsPane.params)

    def newWorkspace(self):
        result = self.maybeClose()
        if result is not None:
            self.mainLogger.info("CLEARING WORKSPACE")
            self.mdiArea.closeAllSubWindows()
            self.start = time.perf_counter()
            self.guiSet()

    def restore_tree(self, parent, datastream, num_childs=None):
        if not num_childs:
            num_childs = datastream.readUInt32()
        for i in range(0, num_childs):
            child = QStandardItem()
            child.read(datastream)
            parent.appendRow(child)
            num_childs = datastream.readUInt32()
            if num_childs > 0:
                self.restore_tree(child, datastream, num_childs)

    def openWorkspace(self):
        sel = QFileDialog.getOpenFileName(parent=self, caption="Open Workspace", directory=QDir.homePath(),
                                          filter="Linnaeo Workspace (*.lno);;Any (*)")
        if sel != ('', ''):
            filename = sel[0]
            self.mainLogger.debug("Opening file: " + str(filename))
            file = QFile(filename)
            file.open(QIODevice.ReadOnly)
            fstream = QDataStream(file)
            self.mainLogger.info("Starting restore")
            sequences = fstream.readQVariantHash()
            titles = fstream.readQVariantList()
            windex = fstream.readUInt32()
            windows = {}
            self.mdiArea.closeAllSubWindows()
            newBModel = widgets.ItemModel(windows, seqTree=True)
            newPModel = widgets.ItemModel(windows)
            self.restore_tree(newBModel.invisibleRootItem(), fstream)
            self.restore_tree(newPModel.invisibleRootItem(), fstream)
            self.start = time.perf_counter()
            self.guiSet(trees=[newBModel, newPModel], data=[sequences, titles, windex])

    def importSequence(self):
        #diag = QFileDialog(self)
        #diag.setDirectory(QDir.homePath())
        #diag.setFileMode(QFileDialog.ExistingFiles)
        #diag.setNameFilter("Fasta File (*.fa *.fasta);;Any (*)")
        sel = QFileDialog.getOpenFileNames(parent=self, caption="Load a Single Sequence", directory=QDir.homePath(),
                                           filter="Fasta File (*.fa *.fasta)")
        print(sel)
        filenames = sel[0]
        afilter = sel[1][-8:-1]
        stype = None
        if afilter == "*.fasta":
            stype = 'fasta'
        for filename in filenames:
            print(filename)
            try:
                with open(filename, 'r'):
                    seq = SeqIO.read(filename, stype)
                    print(seq)
                    seqr = models.SeqR(seq.seq, name=seq.name, id=seq.id, description=seq.description)
                    print([seqr.seq])
                    if [seqr.seq] not in self.sequences.values():
                        self.seqInit(seqr)
                    else:
                        self.mainStatus.showMessage("You've already loaded that sequence!", msecs=4000)
            except:
                if stype == 'fasta':
                    self.mainStatus.showMessage("ERROR -- Please check file. Could this have been an alignment?", msecs=4000)
                else:
                    self.mainStatus.showMessage("ERROR -- Please check file", msecs=3000)

    def importAlignment(self):
        sel = QFileDialog.getOpenFileNames(parent=self, caption="Load an Alignment", directory=QDir.homePath(),
                                           filter="Clustal File (*.aln *.clustal);;Fasta File (*.fa *.fasta);;Any (*)")
        filenames = sel[0]
        afilter = sel[1][-8:-1]
        atype = None
        if afilter == "*.fasta":
            atype = 'fasta'
        elif afilter == "clustal":
            atype = 'clustal'
        count = 0
        for filename in filenames:
            try:
                with open(filename, 'r'):
                    aln = AlignIO.read(filename, atype)
                items = {}
                combo = []

                for seqr in aln:
                    # Go through each sequence in the alignment and import it if it does not already exist.
                    cleanseq = str(seqr.seq).replace('-', '')
                    cleanseq2 = Seq(cleanseq, generic_protein)
                    seqr_clean = models.SeqR(cleanseq2)
                    seqr_clean.convert(seqr)
                    items[seqr_clean.name] = str(seqr.seq)
                    combo.append(seqr_clean)
                    if seqr_clean not in self.sequences.values():
                        num = "" if count == 0 else " "+str(count+1)
                        self.seqInit(seqr_clean, folder='Import%s' % num)
                combo.sort()
                if combo not in self.sequences.values():
                    wid = str(int(self.windex) + 1)
                    self.sequences[wid] = combo
                    sub = self.makeNewWindow(wid, items)
                    self.openWindow(sub)
                    self.windex = self.windex + 1

                else:
                    self.mainStatus.showMessage("Alignment already imported!", msecs=3000)
            except:
                traceback.print_exc()
                self.mainStatus.showMessage("ERROR -- Please check file", msecs=3000)
            finally:
                count += 1

    def exportSequence(self):
        index = self.bioTree.selectedIndexes()
        if index:
            sel = QFileDialog.getSaveFileName(parent=self, caption="Save Alignment", directory=QDir.homePath(),
                                              filter="Fasta File (*.fa *.fasta);;Any (*)")
            if sel != ('', ''):
                afilter = sel[1][-8:-1]
                atype = None
                if afilter == "*.fasta":
                    atype = 'fasta'
                filename = sel[0]# + '.' + atype
                seqR = self.bioModel.itemFromIndex(index[0]).data(role=self.SequenceRole)[0]
                print(seqR)
                with open(filename, 'w'):
                    SeqIO.write(seqR, filename, atype)
                    print(seqR)
        else:
            self.mainStatus.showMessage("Hey, pick a sequence first!", msecs=4000)

    def exportAlignment(self):
        index = self.projectTree.selectedIndexes()
        if index:
            sel = QFileDialog.getSaveFileName(parent=self, caption="Save Alignment", directory=QDir.homePath(),
                                              filter="Clustal File (*.aln *.clustal);;Fasta File (*.fa *.fasta)")
            if sel != ('', ''):
                afilter = sel[1][-8:-1]
                atype = None
                if afilter == "*.fasta":
                    atype = 'fasta'
                elif afilter == "clustal":
                    atype = 'clustal'
                filename = sel[0]# + '.' + atype if sel[0][-1 * len(atype):] != atype else sel[0]
                aliR = self.projectModel.itemFromIndex(index[0]).data(role=self.SequenceRole)
                with open(filename, 'w'):
                    AlignIO.write(aliR, filename, atype)
                    print(aliR)
        else:
            self.mainStatus.showMessage("Hey, pick an alignment first!", msecs=4000)

    def saveWorkspace(self):
        sel = QFileDialog.getSaveFileName(parent=self, caption="Save Workspace", directory=QDir.homePath(),
                                          filter="Linnaeo Workspace (*.lno);;Any (*)")
        filename = sel[0]
        if filename:
            if filename[-4:] != ".lno":
                filename = str(filename + ".lno")
            file = QFile(filename)
            file.open(QIODevice.WriteOnly)
            out = QDataStream(file)
            self.mainLogger.debug("Beginning file save")
            out.writeQVariantHash(self.sequences)
            out.writeQVariantList(self.titles)
            out.writeUInt32(self.windex)
            #####################################
            # Sequence View
            ###
            self.mainLogger.debug("SEQUENCE TREE")
            self.mainLogger.debug("Children count: %s" % self.bioModel.invisibleRootItem().rowCount())
            out.writeUInt32(self.bioModel.invisibleRootItem().rowCount())
            for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
                self.mainLogger.debug("Node: %s" % str(node.text()))
                node.write(out)
                out.writeUInt32(node.rowCount())
            #####################################
            # Alignment View
            # Does not save any metadata! Only the two sequences
            # So I can't save window options at the moment.
            # TODO: Consider adding "window modes" to the node.
            self.mainLogger.debug("ALIGNMENT TREE")
            self.mainLogger.debug("Children count: %s" % self.bioModel.invisibleRootItem().rowCount())
            out.writeUInt32(self.projectModel.invisibleRootItem().rowCount())
            for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
                self.mainLogger.debug("Node: %s" % str(node.text()))
                node.write(out)
                out.writeUInt32(node.rowCount())
            self.mainLogger.debug("Save complete")
            file.flush()
            file.close()
            return True
        else:
            self.mainLogger.debug("No filename chosen; canceling")
            return False

    def get_UniprotId(self):
        if self._currentWindow:
            found = False
            wid = self._currentWindow.wid
            for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
                if node.data(self.WindowRole) == wid:
                    print("Calling UNIPROT search")
                    worker = utilities.GetPDBThread([node.data(role=self.SequenceRole), [node.index()]], parent=self)
                    worker.finished.connect(self.pdbSearchDone)
                    worker.start()
                    found = True
            if not found:
                for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
                    if node.data(self.WindowRole) == wid:
                        print("Calling UNIPROT search")
                        seqs = []
                        # ONLY RUN ON THE TOP SEQUENCE FOR NOW
                        for x in node.data(role=self.SequenceRole):
                            seqs.append(x)
                        worker = utilities.GetPDBThread([[seqs[0]], [node.index()]],
                                                        parent=self)
                        worker.finished.connect(self.pdbSearchDone)
                        worker.start()

        #if self.lastClickedTree == self.bioTree:
        #    indices, seqs = utilities.nodeSelector(self.bioTree, self.bioModel)


            #worker = utilities.GetPDBThread([seqs, indices], parent=self)
            #worker.finished.connect(self.pdbSearchDone)
            #worker.start()


    def pdbSearchDone(self, result):
        if result:
            self.mainLogger.info('Successfully collected PDB from SwissModel repository')
            for match in result:
                worker = utilities.DSSPThread(match, parent=self)
                worker.finished.connect(self.ssDone)
                self.mainLogger.debug('Fowarding structure to DSSP')
                worker.start()
        else:
            self.mainStatus.showMessage('Sorry, no models found for that query! Please manually edit the \
            search by right clicking the node! ', msecs=5000)
        
    def ssDone(self, result):
        print("DSSP FINISHED")
        if result[2]:
            node = self.bioModel.itemFromIndex(result[2])
            if not node:
                node = self.projectModel.itemFromIndex(result[2])
            node.setData(result[0], role=self.StructureRole)
            if node.data(role=self.WindowRole):
                try:
                    sub = self.windows[node.data(role=self.WindowRole)]
                    sub.addStructure(result[0], result[1])
                    if sub == self._currentWindow:

                        self.optionsPane.structureActivate(True)
                except KeyError:
                    self.mainStatus.showMessage("Please open the sequence first!", msecs=3000)


    def copyOut(self):
        if self.lastClickedTree == self.bioTree:
            indices, seqs = utilities.nodeSelector(self.bioTree, self.bioModel)
            fseqs = []
            for seqr in seqs:
                print(seqr.format("fasta"))
                fseqs.append(seqr.format("fasta"))
            QApplication.clipboard().setText("".join(fseqs))
            if len(seqs) > 1:
                self.mainStatus.showMessage("Copied sequences to clipboard!", msecs=1000)
            elif len(seqs) == 1:
                self.mainStatus.showMessage("Copied sequence to clipboard!", msecs=1000)

        elif self.lastClickedTree == self.projectTree:
            self.mainStatus.showMessage("Please use File>Export to save alignments!", msecs=1000)

    def pasteInto(self):
        """
        Only accepts single FASTA paste right now
        """
        # FASTA DETECTION
        if self.lastClickedTree == self.bioTree:
            nline = None
            sid = None
            desc = None
            name = None
            bars = []
            spaces = []
            try:
                clip = QApplication.clipboard().text()
                if clip[0] == ">":
                    for index in range(len(clip)):
                        if not nline and clip[index] == "\n":
                            nline = clip[:index]
                            seq = clip[index + 1:]
                    seq = seq.replace('\n', '')
                    seq = seq.replace('\r', '')
                    bseq = MutableSeq(seq)

                    for index in range(len(nline)):
                        if nline[index] == " ":
                            spaces.append(index)
                        if nline[index] == "|":
                            bars.append(index)
                    if len(bars) > 0:
                        sid = nline[bars[0] + 1:bars[1]]
                        self.mainLogger.debug("Extracted ID for sequence: %s" % sid)
                        desc = nline[bars[1] + 1:]
                    elif not sid:
                        sid = nline[1:9]
                        name = sid

                    if sid and bseq:
                        seqr = models.SeqR(bseq, id=sid, name=sid)
                        if desc:
                            seqr.description = desc
                        self.seqInit(seqr)
                else:
                    self.mainStatus.showMessage("Please only paste in FASTA format!", msecs=1000)
            except:
                self.mainStatus.showMessage("Please only paste in FASTA format!", msecs=1000)
        else:
            self.mainStatus.showMessage("Please choose paste destination", msecs=1000)

    def addFolder(self):
        new = QStandardItem("New Folder")
        new.setData("New Folder")
        try:
            node = self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0])
            if not node.data(role=self.WindowRole):
                node.appendRow(new)
            else:
                self.bioModel.appendRow(new)
        except IndexError:
            try:
                node = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
                if not node.data(role=self.WindowRole):
                    node.appendRow(new)
                else:
                    self.projectModel.appendRow(new)
            except IndexError:
                if self.lastClickedTree:
                    self.mainLogger.debug("Nothing selected so adding folder at last click")
                    self.lastClickedTree.model().appendRow(new)
                else:
                    self.mainStatus.showMessage("Select a location first!", msecs=1000)

    def deleteNode(self):
        for index in self.lastClickedTree.selectedIndexes():
            node = self.lastClickedTree.model().itemFromIndex(index)
            wid = node.data(role=self.WindowRole)
            try:
                sub = self.windows[wid]
                if sub:
                    self.mainLogger.debug("Deleting node from tree: " + str(sub.windowTitle()))
                    sub.close()
                    self.windows.pop(wid)
                    self.sequences.pop(wid)
                    self.pruneNames()
                    self.bioModel.updateNames(self.titles)
                    self.bioModel.updateWindows(self.windows)
                    self.projectModel.updateNames(self.titles)
                    self.projectModel.updateWindows(self.windows)
            except KeyError:
                pass
            except ValueError:
                pass
            self.lastClickedTree.model().removeRow(index.row(), index.parent())
            if self.mdiArea.tabbed:
                self.mdiArea.activeSubWindow().showMaximized()

    def tileWindows(self):
        if self.mdiArea.tabbed:
            self.mdiArea.toggleTabs()
        self.mdiArea.tileSubWindows()

    def cascadeWindows(self):
        if self.mdiArea.tabbed:
            self.mdiArea.toggleTabs()
        self.mdiArea.cascadeSubWindows()

    def closeTab(self):
        """ This duplicates inherent functionality of QMDISubWindow, so I can't add shortcut to menu"""
        try:
            self.mdiArea.activeSubWindow().close()
        except AttributeError:
            pass

    def closeAllTabs(self):
        self.mdiArea.closeAllSubWindows()

    def dupeNameMsg(self):
        self.mainStatus.showMessage("Please choose a unique name!", msecs=1000)


class Debug:
    def DEBUG(self):
        # STRICTLY FOR TESTING -- FAKE DATA.
        test1 = ['GPI1A', 'MSLSQDATFVELKRHVEANEKDAQLLELFEKDPARFEKFTRLFATPDGDFLFDF' +
                 'SKNRITDESFQLLMRLAKSRGVEESRNAMFSAEKINFTENRAVLHVALRNRANRP' +
                 'ILVDGKDVMPDVNRVLAHMKEFCNEIISGSWTGYTGKKITDVVNIGIGGSDLGPL' +
                 'MVTESLKNYQIGPNVHFVSNVDGTHVAEVTKKLNAETTLFIIASKTFTTQETITN' +
                 'AETAKEWFLAKAGDAGAVAKHFVALSTNVTKAVEFGIDEKNMFEFWDWVGGRYSL' +
                 'WSAIGLSIAVHIGFDNYEKLLDGAFSVDEHFVNTPLEKNIPVILAMIGVLYNNIY' +
                 'GAETHALLPYDQYMHRFAAYFQQGDMESNGKFVTRHGQRVDYSTGPIVWGEPGTN' +
                 'GQHAFYQLIHQGTRLIPADFIAPVKTLNPIRGGLHHQILLANFLAQTEALMKGKT' +
                 'AAVAEAELKSSGMSPESIAKILPHKVFEGNKPTTSIVLPVVTPFTLGALIAFYEH' +
                 'KIFVQGIIWDICSYDQWGVELGKQLAKVIQPELASADTVTSHDASTNGLIAFIKNNA']
        seq_GPI1A = MutableSeq(test1[1], generic_protein)
        gpi1a = models.SeqR(seq_GPI1A, id='Q9U1Q2', name=test1[0])
        #gpi1a.id = test1[0]
        test2 = ['GPI1B', 'MIFELFRFIFRKKKMLGYLSDLIGTLFIGDSTEKAMSLSQDATFVELKRHVEANE' +
                 'KDAQLLELFEKDPARFEKFTRLFATPDGDFLFDFSKNRITDESFQLLMRLAKSRG' +
                 'VEESRNAMFSAEKINFTENRAVLHVALRNRANRPILVDGKDVMPDVNRVLAHMKE' +
                 'FCNEIISGSWTGYTGKKITDVVNIGIGGSDLGPLMVTESLKNYQIGPNVHFVSNV' +
                 'DGTHVAEVTKKLNAETTLFIIASKTFTTQETITNAETAKEWFLAKAGDAGAVAKH' +
                 'FVALSTNVTKAVEFGIDEKNMFEFWDWVGGRYSLWSAIGLSIAVHIGFDNYEKLL' +
                 'DGAFSVDEHFVNTPLEKNIPVILAMIGVLYNNIYGAETHALLPYDQYMHRFAAYF' +
                 'QQGDMESNGKFVTRHGQRVDYSTGPIVWGEPGTNGQHAFYQLIHQGTRLIPADFI' +
                 'APVKTLNPIRGGLHHQILLANFLAQTEALMKGKTAAVAEAELKSSGMSPESIAKI' +
                 'LPHKVFEGNKPTTSIVLPVVTPFTLGALIAFYEHKIFVQGIIWDICSYDQWGVEL' +
                 'GKQLAKVIQPELASADTVTSHDASTNGLIAFIKNNA']
        seq_GPI1B = MutableSeq(test2[1], generic_protein)
        gpi1b = models.SeqR(seq_GPI1B, id='Q7K707', name=test2[0])
        test = [gpi1a, gpi1b]
        for i in test:
            self.seqInit(i)

    def _queryTrees(self):
        """ Currently not used. Spits out everything in the tree.
        """
        print("\n\nBIOROOT")
        for child in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            print("Text: ", child.text())
            print("Name: ", child.data(role=Qt.UserRole + 1))
            print("Seq: ", child.data(role=Qt.UserRole + 2))
            print("Window Index: ", child.data(role=Qt.UserRole + 3))
        print("ALIGNMENT ROOT")
        for child in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
            print("Text: ", child.text())
            print("Name: ", child.data(role=Qt.UserRole + 1))
            print("Seq: ", child.data(role=Qt.UserRole + 2))
            print("Window Index: ", child.data(role=Qt.UserRole + 3))

        """
        # SAVED THIS FOR LATER!
        # Read in config (linux)
        config = configparser.ConfigParser()
        try:
            config.read(str(Path.home())+"/.linnaeo/config.ini")
            # Open last used workspace automatically.
            if config['RECENTS']['LAST'] != "":
                last = config['RECENTS']['LAST']
                f = QFile(last)
                f.open(QIODevice.ReadOnly)
                model = workspace.WorkspaceModel()
                f.close()
                self.workspaceTree.setModel(model)
        except:
            print("No config file found!")
            """
