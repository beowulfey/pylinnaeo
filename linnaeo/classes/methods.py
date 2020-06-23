import time
import traceback

from Bio import SeqIO, AlignIO
from Bio.Alphabet import generic_protein
from Bio.Seq import MutableSeq, Seq
from PyQt5.QtCore import QFile, QIODevice, QDataStream, Qt, QDir, QUrl
from PyQt5.QtGui import QStandardItem, QFontMetricsF, QIcon, QDesktopServices
from PyQt5.QtWidgets import QFileDialog, QApplication, qApp

from linnaeo.classes import widgets, models, utilities
from linnaeo.classes.displays import AboutDialog
from linnaeo.classes.utilities import lookupTheme


class Slots:
    """
    Split out methods for clarity. These are all simple slot functions associated with UI actions!
    """

    #def increaseTextSize(self):
    #    self._currentWindow.widget().increaseFont()
    #    self.mainStatus.showMessage("Setting font size to %s " % self._currentWindow.widget().getFontSize(), msecs=1000)

    #def decreaseTextSize(self):
    #    self._currentWindow.widget().decreaseFont()
    #    self.mainStatus.showMessage("Setting font size to %s " % self._currentWindow.widget().getFontSize(), msecs=1000)

    def setCurrentWindow(self):
        """
        Important function for changing the currently displayed window. It updates the stored "current window",
        and modifies the options pane depending on the status of the new window (such as disabling buttons, resetting
        the reference sequence, etc
        """
        self._currentWindow = self.mdiArea.activeSubWindow()
        if self._currentWindow:
            refseq = self._currentWindow.widget().refseq
            self.optionsPane.comboReference.clear()
            self.optionsPane.comboReference.addItem("Select seq...")
            self.optionsPane.comboReference.setCurrentIndex(0)
            self._currentWindow.widget().refseq = refseq
            del refseq
            for name in self._currentWindow.widget().splitNames:
                self.optionsPane.comboReference.addItem(name)
            del name
            if self._currentWindow.widget().refseq is not None:
                for x in range(self.optionsPane.comboReference.count()):
                    if self.optionsPane.comboReference.itemText(x) == self._currentWindow.widget().splitNames[
                            self._currentWindow.widget().refseq]:
                        self.optionsPane.comboReference.setCurrentIndex(x)
                    del x
            if not self.optionsPane.checkStructure.isEnabled():
                if self._currentWindow.widget().dssps:
                    self.optionsPane.structureActivate(True)
                else:
                    if not self.optionsPane.buttonStructure.isEnabled() and self._currentWindow not in self.runningDSSP:
                        self.optionsPane.buttonStructure.setEnabled(True)
            elif not self._currentWindow.widget().dssps:
                self.optionsPane.structureActivate(False)
        else:
            self.optionsPane.structureActivate(False)
            self.optionsPane.buttonStructure.setEnabled(False)

    def showAbout(self):
        """ Open ABOUT dialog window. """
        qDialog = AboutDialog(self)
        qDialog.exec()

    def saveImage(self):
        """ Captures the current alignment window in entirety. Very crude. """
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
            del file, w, sw

    def toggleOptionsPane(self, state):
        """ Permanent button for hiding or displaying the right-side options pane for alignment controls. """
        self.optionsPane.show() if state else self.optionsPane.hide()
        del state

    def toggleRuler(self, state):
        """ Turn horizontal ruler on/off """
        if self._currentWindow:
            self._currentWindow.widget().toggleRuler(state)
        del state

    def toggleColors(self, state):
        """ Turn sequence colors on/off"""
        if self._currentWindow:
            self._currentWindow.widget().toggleColors(state)
        del state

    def toggleStructure(self, state):
        """ Turn structure visible on/off """
        if self._currentWindow:
            self._currentWindow.widget().toggleStructure(bool(state))
        del state

    def toggleConsv(self, state):
        """ Turn reference sequence on/off """
        if self._currentWindow:
            self._currentWindow.widget().setConsvColors(state)
        del state

    def selectReference(self, index):
        """ Simple forwarding of selected reference sequence for coloring of sequence. """
        self._currentWindow.widget().setReference(self.optionsPane.comboReference.itemText(index))
        del index

    def changeTheme(self):
        """ Simple forwarding of selected theme. """
        self.colorPane.clear()
        desc = lookupTheme(self.optionsPane.comboTheme.currentText()).getDesc()
        self.colorPane.insertHtml(desc)
        font = self.colorPane.document().defaultFont()
        fmF = QFontMetricsF(font)
        text = self.colorPane.document().toPlainText()
        textSize = fmF.size(0, text)
        hgt = textSize.height()+4
        self.colorPane.setFixedHeight(hgt)
        del desc, font, fmF, text, textSize, hgt
        if self.optionsPane.comboTheme.currentText() == "Conservation":
            self.optionsPane.checkConsv.setChecked(True)
        if self._currentWindow:
            self.optionsPane.comboReference.setCurrentIndex(1)
            self._currentWindow.widget().setTheme(self.optionsPane.comboTheme.currentText())


    def changeFont(self, font):
        """
        Forwards selected font to the widget. If the chosen font is too different from the width of my symbol font,
        it declares it incompatible and disables structure visualization.
        """
        if self._currentWindow:
            self._currentWindow.widget().setFont(font)
            fmF = QFontMetricsF(font)
            if not fmF.averageCharWidth() - 0.5 <= self._currentWindow.widget().ssFontWidth <= \
                   fmF.averageCharWidth() + 0.5:
                self.mainLogger.debug("Font size is too different from symbol font; disabling structure!")
                self.mainStatus.showMessage("Font incompatible with structure symbols; please choose another!",
                                            msecs=6000)
                self.optionsPane.structureActivate(False)
            elif self._currentWindow.widget().dssps and not self.optionsPane.checkStructure.isEnabled():
                self.optionsPane.structureActivate(True)
            del fmF
        del font

    def changeFontSize(self, size):
        """ Simple forwarding of font size to the widget. """
        if self._currentWindow:
            self._currentWindow.widget().setFontSize(size)
        del size

    def refreshParams(self, window):
        """ Simple; push the current window options from the options pane to a window on selection. """
        if window:
            window.widget().setParams(self.optionsPane.params)
            del window

    def showColorDesc(self, state):
        """if state:
            self.gridLayout_2.addWidget(self.colorPane, 1, 0)
            self.colorPane.show()
        else:
            self.gridLayout_2.removeWidget(self.colorPane)
            self.colorPane.hide()"""
        if state:
            self.optionsPane.verticalLayout.insertWidget(12,self.colorPane)
            self.colorPane.show()
        else:
            self.optionsPane.verticalLayout.removeWidget(self.colorPane)
            self.colorPane.hide()

    def openThemeHelp(self):
        QDesktopServices.openUrl(QUrl('https://beowulfey.github.io/linnaeo/linnaeo/resources/docs/themes'))

    def restore_tree(self, parent, datastream, num_childs=None):
        """ Function that acts during opening a workspace. Rebuilds a tree by iterating through it. """
        if not num_childs:
            num_childs = datastream.readUInt32()
        for i in range(0, num_childs):
            child = QStandardItem()
            child.read(datastream)
            parent.appendRow(child)
            num_childs = datastream.readUInt32()
            if num_childs > 0:
                self.restore_tree(child, datastream, num_childs)
            del child
        del parent, datastream, num_childs

    def rebuildTrees(self):
        """
        Run after loading a saved file.
        At this point, the sequences and the alignments both have Window IDs applied -- but the windows
        no longer exist. This rebuilds the windows using the saved window IDs and updates the respective models.
        """
        for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
            if node.data(role=self.SequenceRole):
                self.mainLogger.info("Loading sequence: "+node.data())
                self.mainStatus.showMessage("Loading sequence: %s" % node.data(), msecs=3000)
                ali = {}  # empty dict needed to send to open window
                wid = node.data(role=self.WindowRole)
                seqr = node.data(role=self.SequenceRole)[0]
                ali[seqr.name] = str(seqr.seq)
                self.makeNewWindow(wid, ali, nonode=True)

        for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
            if node.data(role=self.SequenceRole):
                self.mainLogger.info("Loading alignment: "+node.data())
                self.mainStatus.showMessage("Loading alignment: %s" % node.data(), msecs=3000)
                seqs = {}
                wid = node.data(role=self.WindowRole)
                seqr = node.data(role=self.SequenceRole)
                for seq in seqr:
                    seqs[seq.name] = str(seq.seq)
                worker = utilities.AlignThread(self, seqs, seqtype=3, num_threads=self.threadpool.maxThreadCount())
                worker.start()
                worker.finished.connect(worker.deleteLater)
                worker.finished.connect(worker.quit)
                worker.wait()
                ali = worker.aligned
                del worker
                self.makeNewWindow(wid, ali, nonode=True)
        self.bioModel.updateWindows(self.windows)
        self.projectModel.updateWindows(self.windows)
        self.mainStatus.showMessage(f"Loading complete! took {float(time.perf_counter()-self.start):.2f} seconds", msecs=4000)
        self.mainLogger.debug("Regenerating windows took took %f seconds" % float(time.perf_counter() - self.start))
        del node

    def saveWorkspace(self):
        """
        Saves the current workspace as is. Saves the trees and all of the nodes (and node data). Does not save
        windows or window options, but does save the structure information.
        """
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
            del sel, filename, file, out, node
            return True
        else:
            self.mainLogger.debug("No filename chosen; canceling")
            del sel
            return False

    def importSequence(self):
        """ Uses BioPython to parse an opened file."""
        sel = QFileDialog.getOpenFileNames(parent=self, caption="Load a Single Sequence", directory=QDir.homePath(),
                                           filter="Fasta File (*.fa *.fasta)")
        filenames = sel[0]
        afilter = sel[1][-8:-1]
        stype = None
        if afilter == "*.fasta":
            stype = 'fasta'
        for filename in filenames:
            try:
                with open(filename, 'r'):
                    seq = SeqIO.read(filename, stype)
                    print(seq)
                    self.mainStatus.showMessage("Loading sequence: %s" % seq.id, msecs=3000)
                    seqr = models.SeqR(seq.seq, name=seq.name, id=seq.id, description=seq.description)
                    if [seqr.seq] not in self.sequences.values():
                        self.seqInit(seqr)
                    else:
                        self.mainStatus.showMessage("You've already loaded that sequence!", msecs=4000)
                    del seq, seqr
            except:
                if stype == 'fasta':
                    self.mainStatus.showMessage("ERROR -- Please check file. Could this have been an alignment?",
                                                msecs=4000)
                else:
                    self.mainStatus.showMessage("ERROR -- Please check file", msecs=3000)
                    traceback.print_exc()
        del sel, filenames, stype, afilter

    def importAlignment(self):
        """
        Uses BioPython to parse an alignment file. Imports all the subcomponents as Sequences, and adds
        those to the sequence tree under a folder called 'Import'. If you import two alignments at once, it will add
        to additional folders (Import 2, Import 3, etc).
        """
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
                        num = "" if count == 0 else " " + str(count + 1)
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
                del aln, combo, items, seqr, cleanseq, cleanseq2, seqr_clean, num, wid, sub
            except:
                traceback.print_exc()
                self.mainStatus.showMessage("ERROR -- Please check file", msecs=3000)
            finally:
                count += 1
        del sel, filenames, afilter, count


    def exportSequence(self):
        """ Simple method to export a sequence. Currently only does Fasta but it would be easy to implement more. """
        index = self.bioTree.selectedIndexes()
        if index:
            sel = QFileDialog.getSaveFileName(parent=self, caption="Save Alignment", directory=QDir.homePath(),
                                              filter="Fasta File (*.fa *.fasta);;Any (*)")
            if sel != ('', ''):
                afilter = sel[1][-8:-1]
                atype = None
                if afilter == "*.fasta":
                    atype = 'fasta'
                filename = sel[0]  # + '.' + atype
                seqR = self.bioModel.itemFromIndex(index[0]).data(role=self.SequenceRole)[0]
                #print(seqR)
                with open(filename, 'w'):
                    SeqIO.write(seqR, filename, atype)
                    print(seqR)
            del index, sel, afilter, atype, filename, seqR
        else:
            self.mainStatus.showMessage("Hey, pick a sequence first!", msecs=4000)

    def exportAlignment(self):
        """ Method to export an alignment using BioPython. """
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
                filename = sel[0]  # + '.' + atype if sel[0][-1 * len(atype):] != atype else sel[0]
                aliR = self.projectModel.itemFromIndex(index[0]).data(role=self.SequenceRole)
                with open(filename, 'w'):
                    AlignIO.write(aliR, filename, atype)
                    print(aliR)
            del index, sel, afilter, atype, filename, aliR
        else:
            self.mainStatus.showMessage("Hey, pick an alignment first!", msecs=4000)

    def statuslogger(self, status):
        """ Mini function to show statusbar messages from other threads. """
        self.mainStatus.showMessage(status,msecs=5000)
        del status

    def get_UniprotId(self):
        """
        Uses BioServices to query Uniprot using the saved Sequence ID. See the GetPDBThread for more info.
        For some reason, on macs the button fails to be disabled after clicking when I freeze the binary... no idea why.
        """
        if self._currentWindow:
            self.runningDSSP.append(self._currentWindow)
            self.optionsPane.buttonStructure.setEnabled(False)
            found = False
            wid = self._currentWindow.wid
            for node in utilities.iterTreeView(self.bioModel.invisibleRootItem()):
                if node.data(self.WindowRole) == wid:
                    self.mainStatus.showMessage("Beginning PDB search via UNIPROT")
                    self.mainLogger.info("Running UNIPROT search with the sequence")
                    worker = utilities.GetPDBThread([node.data(role=self.SequenceRole), [node.index()]], parent=self)
                    worker.logger.connect(self.statuslogger)
                    worker.finished.connect(self.pdbSearchDone)
                    worker.finished.connect(worker.deleteLater)
                    worker.finished.connect(worker.quit)
                    worker.start()
                    del worker, node
                    found = True
            if not found:
                for node in utilities.iterTreeView(self.projectModel.invisibleRootItem()):
                    if node.data(self.WindowRole) == wid:
                        self.mainStatus.showMessage("Beginning PDB search via UNIPROT")
                        self.mainLogger.info("Running UNIPROT search with first sequence of the alignment")
                        seqs = []
                        # ONLY RUN ON THE TOP SEQUENCE FOR NOW
                        for x in node.data(role=self.SequenceRole):
                            seqs.append(x)
                        worker = utilities.GetPDBThread([[seqs[0]], [node.index()]],
                                                        parent=self)
                        worker.finished.connect(self.pdbSearchDone)
                        worker.finished.connect(worker.deleteLater)
                        worker.finished.connect(worker.quit)
                        worker.start()
                        del worker, node, x
            del found, wid

    def pdbSearchDone(self, result):
        """ Immediately run after a PDB is acquired from the previous method. Calls DSSP on the saved PDB object. """
        if result:
            self.mainLogger.info('Successfully collected PDB from SwissModel repository')
            self.mainStatus.showMessage('Successfully collected PDB from SwissModel',msecs=5000)
            for match in result:
                worker = utilities.DSSPThread(match, parent=self)
                self.mainLogger.debug('Fowarding structure to DSSP')
                worker.finished.connect(self.ssDone)
                worker.finished.connect(worker.deleteLater)
                worker.finished.connect(worker.quit)
                worker.start()
                del worker
            del match, result
        else:
            self.optionsPane.buttonStructure.setEnabled(True)
            self.mainStatus.showMessage('SORRY--QUERY FAILED! No model was found!', msecs=5000)
            del result

    def ssDone(self, result):
        """ After running DSSP, this adds the DSSP data to the node that corresponds to the window you ran it on. """
        if result[2]:
            self.mainLogger.info("DSSP run completed successfully")
            self.mainStatus.showMessage("DSSP run completed successfully", msecs=5000)
            #print(result[2])
            node = self.bioModel.itemFromIndex(result[2])
            if not node:
                node = self.projectModel.itemFromIndex(result[2])
            node.setData(result[0], role=self.StructureRole)
            if node.data(role=self.WindowRole):
                try:
                    sub = self.windows[node.data(role=self.WindowRole)]
                    self.mainLogger.info("Adding DSSP to subwindow")
                    self.mainStatus.showMessage("DSSP found; adding to window %s" % node.text(), msecs=4000)
                    sub.addStructure(result[0], result[1])
                    if sub in self.runningDSSP:
                        self.runningDSSP.remove(sub)
                    if sub == self._currentWindow:
                        self.optionsPane.structureActivate(True)
                    del sub
                except KeyError:
                    self.mainStatus.showMessage("Please open the sequence first!", msecs=3000)
            del node
        else:
              self.mainLogger.info("DSSP failed -- please check your sequence ID and that DSSP is available")
              self.mainStatus.showMessage("DSSP failed -- please check your sequence ID and that DSSP is available",msecs=5000)
        del result

    def copyOut(self):
        """ Copies the sequence from the selected sequence in the Biotree. Does not copy alignments. """
        if self.lastClickedTree == self.bioTree:
            indices, seqs = utilities.nodeSelector(self.bioTree, self.bioModel)
            fseqs = []
            for seqr in seqs:
                print(seqr.format("fasta"))
                fseqs.append(seqr.format("fasta"))
            QApplication.clipboard().setText("".join(fseqs))
            if len(seqs) > 1:
                self.mainStatus.showMessage("Copied sequences to clipboard!", msecs=2000)
            elif len(seqs) == 1:
                self.mainStatus.showMessage("Copied sequence to clipboard!", msecs=2000)
            del indices, seqs, fseqs
        elif self.lastClickedTree == self.projectTree:
            self.mainStatus.showMessage("Please use File>Export to save alignments!", msecs=4000)

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
                        del seqr
                    del index, seq, bseq
                else:
                    self.mainStatus.showMessage("Please only paste in FASTA format!", msecs=2000)
                del clip
            except:
                self.mainStatus.showMessage("Please only paste in FASTA format!", msecs=2000)
            del nline, sid, desc, name, bars, spaces
        else:
            self.mainStatus.showMessage("Please choose paste destination", msecs=2000)

    def addFolder(self):
        """ Adds a folder to the currently selected tree. """
        new = QStandardItem("New Folder")
        new.setData("New Folder")
        try:
            node = self.bioModel.itemFromIndex(self.bioTree.selectedIndexes()[0])
            if not node.data(role=self.WindowRole):
                node.appendRow(new)
            else:
                self.bioModel.appendRow(new)
            del node
        except IndexError:
            try:
                node = self.projectModel.itemFromIndex(self.projectTree.selectedIndexes()[0])
                if not node.data(role=self.WindowRole):
                    node.appendRow(new)
                else:
                    self.projectModel.appendRow(new)
                del node
            except IndexError:
                if self.lastClickedTree:
                    self.mainLogger.debug("Nothing selected so adding folder at last click")
                    self.lastClickedTree.model().appendRow(new)
                else:
                    self.mainStatus.showMessage("Select a location first!", msecs=1000)
        del new

    def deleteNode(self):
        for index in self.lastClickedTree.selectedIndexes():
            node = self.lastClickedTree.model().itemFromIndex(index)
            wid = node.data(role=self.WindowRole)
            try:
                sub = self.windows[wid]
                if sub:
                    self.mainLogger.debug("Deleting node from tree: " + str(sub.windowTitle()))
                    sub.delete()
                    self.windows.pop(wid)
                    self.sequences.pop(wid)
                    self.pruneNames()
                    self.bioModel.updateNames(self.titles)
                    self.bioModel.updateWindows(self.windows)
                    self.projectModel.updateNames(self.titles)
                    self.projectModel.updateWindows(self.windows)
                    del sub
            except KeyError:
                pass
            except ValueError:
                pass
            self.lastClickedTree.model().removeRow(index.row(), index.parent())
            if self.mdiArea.tabbed and self._currentWindow:
                self.mdiArea.activeSubWindow().showMaximized()
            del index, node, wid

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
        # STRICTLY FOR TESTING -- some INIT DATA.
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
        # gpi1a.id = test1[0]
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
        del test1, test2, seq_GPI1A, gpi1a, seq_GPI1B, gpi1b, test, i

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


