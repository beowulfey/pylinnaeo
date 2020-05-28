import time

import Bio
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Seq import MutableSeq
from PyQt5.QtCore import QFile, QIODevice, QDataStream, Qt, QDir
from PyQt5.QtGui import QStandardItem
from PyQt5.QtWidgets import QFileDialog, QApplication

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

    def showAbout(self):
        qDialog = AboutDialog(self)
        qDialog.exec()

    def setSizing(self):
        print("being clicked?", self.beingClicked)
        if not self.beingClicked:
            self.beingClicked = True
            if self._currentWindow and self._currentWindow.isMaximized():  # and self.mdiArea.activeSubWindow().isMaximized():
                print("REDRAWING FRAME FROM MAIN")
                self._currentWindow.widget().userIsResizing = True
                self._currentWindow.widget().seqArrange(color=False) #, rulers=False)
        elif self.beingClicked:
            self.beingClicked = False
            if self._currentWindow and self._currentWindow.isMaximized():
                print("DONE REDRAWING FROM MAIN")
                self._currentWindow.widget().userIsResizing = False
                self._currentWindow.widget().seqArrange()

    def saveImage(self):
        file = QFileDialog.getSaveFileName(self, "Save as...", "name",
                                            "PNG (*.png);; BMP (*.bmp);;TIFF (*.tiff *.tif);; JPEG (*.jpg *.jpeg)");
        self._currentWindow.widget().seqArrange(color=True, rulers=True)
        self._currentWindow.widget().grab().save(file[0]+file[1][-5:-1])

    def toggleRulers(self):
        self._currentWindow.widget().toggleRulers()

    def toggleColors(self):
        self._currentWindow.widget().toggleColors()

    # SINGLE-USE SLOTS
    def newWorkspace(self):
        result = self.maybeClose()
        if result is not None:
            self.mainLogger.info("CLEARING WORKSPACE")
            self.mdiArea.closeAllSubWindows()
            self.start = time.perf_counter()
            self.guiSet()

    def restore_tree(self, parent, datastream, num_childs=None):
        if not num_childs:
            #print("First line: reading UInt32")
            num_childs = datastream.readUInt32()
            #print(num_childs)
        for i in range(0, num_childs):
            #print("Reading node")
            child = QStandardItem()
            child.read(datastream)
            parent.appendRow(child)
            #print(child.data(role=Qt.UserRole + 1))
            #print(child.data(role=Qt.UserRole + 2))
            #print(child.data(role=Qt.UserRole + 3))
            num_childs = datastream.readUInt32()
            if num_childs > 0:
                #print("reading children")
                self.restore_tree(child, datastream, num_childs)

    def openWorkspace(self):
        sel = QFileDialog.getOpenFileName(parent=self, caption="Open Workspace", directory=QDir.homePath(),
                                          filter="Linnaeo Workspace (*.lno);;Any (*)")
        filename = sel[0]
        self.mainLogger.debug("Opening file: " + str(filename))
        file = QFile(filename)
        file.open(QIODevice.ReadOnly)
        fstream = QDataStream(file)
        self.mainLogger.info("Starting restore")
        sequences = fstream.readQVariantHash()
        #print("Sequences: ", sequences)
        titles = fstream.readQVariantList()
        #print("Titles: ", titles)
        windex = fstream.readUInt32()
        #print("Windex: ", windex)
        windows = {}
        self.mdiArea.closeAllSubWindows()
        newBModel = widgets.ItemModel(windows, seqTree=True)
        newPModel = widgets.ItemModel(windows)
        self.restore_tree(newBModel.invisibleRootItem(), fstream)
        self.restore_tree(newPModel.invisibleRootItem(), fstream)
        self.start = time.perf_counter()
        self.guiSet(trees=[newBModel, newPModel], data=[sequences, titles, windex])

    def newSequence(self):
        sel = QFileDialog.getOpenFileName(parent=self, caption="Load a Single Sequence", directory=QDir.homePath(),
                                          filter="Fasta File (*.fa *.fasta);;Any (*)")
        seqr = None
        filename = sel[0]
        badfile = True
        if filename[-3:] == ".fa":
            badfile = False
        elif filename[-6:] == ".fasta":
            badfile = False
        if not badfile:
            try:
                seq = SeqIO.read(filename, "fasta")
                if seq:
                    #print("name: ", seq.name)
                    #print("id: ", seq.id)
                    #print("desc: ", seq.description)
                    seqr = models.SeqR(seq.seq, name=seq.name, id=seq.id, description=seq.description)
                    self.seqInit(seqr)
            except:
                self.mainStatus.showMessage("ERROR -- Please check file", msecs=3000)
        else:
            self.mainStatus.showMessage("Please only add fasta file!", msecs=1000)

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
            #print("TREE DATA")
            #self.queryTrees()
            #print("Writing Sequences, Titles, Windex")
            #print(self.sequences)
            #print(self.titles)
            #print(type(self.windex), self.windex)
            out.writeQVariantHash(self.sequences)
            out.writeQVariantList(self.titles)
            out.writeUInt32(self.windex)
            #####################################
            # Sequence View
            ###
            print("SEQUENCE TREE")
            print("Invs.Root Children: ", self.bioModel.invisibleRootItem().rowCount())
            out.writeUInt32(self.bioModel.invisibleRootItem().rowCount())
            for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
                print("Text: ", str(node.text()))
                print("Data1: ", str(node.data()))
                print("Data2: ", str(node.data(role=self.SequenceRole)))
                print("Data3: ", str(node.data(role=self.WindowRole)))
                node.write(out)
                out.writeUInt32(node.rowCount())
            #####################################
            # Alignment View
            # Does not save any metadata! Only the two sequences
            # So I can't save window options at the moment.
            # TODO: Consider adding "window modes" to the node.
            print("ALIGNMENT TREE")
            print("Invs.Root Children: ", self.bioModel.invisibleRootItem().rowCount())
            out.writeUInt32(self.projectModel.invisibleRootItem().rowCount())
            for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
                print("Text: ", str(node.text()))
                print("Data1: ", str(node.data()))
                print("Data2: ", str(node.data(role=self.SequenceRole)))
                print("Data3: ", str(node.data(role=self.WindowRole)))
                node.write(out)
                out.writeUInt32(node.rowCount())
            self.mainLogger.debug("Save complete")
            file.flush()
            file.close()
            return True
        else:
            #print("No filename chosen; canceling")
            return False

    def copyOut(self):
        if self.lastClickedTree == self.bioTree:
            seqs = []
            indices = self.lastClickedTree.selectedIndexes()
            for index in indices:
                node = self.bioModel.itemFromIndex(index)
                seqr = node.data(role=self.SequenceRole)[0]
                print(seqr.format("fasta"))
                seqs.append(seqr.format("fasta"))
            QApplication.clipboard().setText("".join(seqs))
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
                        self.mainLogger.debug("Extracted ID for sequence: ", sid)
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
        gpi1a = models.SeqR(seq_GPI1A, id=test1[0], name=test1[0])
        gpi1a.id = test1[0]
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
        gpi1b = models.SeqR(seq_GPI1B, id=test2[0], name=test2[0])
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