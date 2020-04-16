#!/usr/bin/python3

# BIOGLOT
# This is my shitty code for managing all the parts I'm going to somehow
# hack together into a piece of software, maybe.

# Included here are my class objects for managing sequences.

class DnaContainer(object):
# Container for any and all DNA sequences. 

    def __init__(self):
        self.sequence = []
        self.readframe = 1
        self.translated_long = []
        self.translated = []
        return

    #def __cleanSequence( ?
    #def __cleanReadframe( ?

    # DEM GET METHODS DOE
    def getSequence(self):
        return self.sequence
    def getReadframe(self):
        return self.readframe
    def getLength(self):
        self.seqlength = len(self.sequence)
        return self.seqlength
    def getTranslationLong(self):
        # Easy -- originally designed codon reader for three letter acronyms. However,
        # I guess I thought the one letter would be easier to read so it's the default. 
        # Soooooo we're ROLLING WITH IT!
        if self.translated_long == []:
            self.__translateSeq()
        return self.translated_long
    
    # ONCE YOU SET YOU WON'T FORGET!
    def setSequence(self, user_sequence):
        self.sequence = user_sequence
        return self.sequence

    def getTranslation(self):
        if self.translated == []:
            if self.translated_long == []:
                self.__translateSeq()
            translong_list = self.translated_long.split()
        
            aa_form_swapper = {
                    'A' : 'ala', 'C' : 'cys', 'D' : 'asp', 'E' : 'glu',
                    'F' : 'phe', 'G' : 'gly', 'H' : 'his', 'I' : 'ile',
                    'K' : 'lys', 'L' : 'leu', 'M' : 'met', 'N' : 'asn',
                    'P' : 'pro', 'Q' : 'gln', 'R' : 'arg', 'S' : 'ser',
                    'T' : 'thr', 'V' : 'val', 'W' : 'trp', 'Y' : 'tyr',
                    '*' : ' * '
                     }

            # Swap out three-letter  for one-letter and join into string.
            translation_list = []
            for long_resi in translong_list:
                for sform, lform in aa_form_swapper.items():
                    if long_resi.lower() in lform:
                        translation_list.append(sform)
            self.translated = ''.join(translation_list)
        return self.translated

    def __translateSeq(self):
        # Not accessible outside: done as part of getTranslation's
        # FIRST PART: make a list of codons based on reading frame
        codon_split = []
        codon_count = int(self.getLength()/3)
        rf_index = self.readframe - 1
        
        for x in range(rf_index,codon_count):
            codon_split.append(self.sequence[rf_index:rf_index + 3])
            rf_index += 3

        # SECOND PART: look up translation
        codon_table = {
                'Ala' : ['GCG','GCC','GCA','GCT'],
                'Cys' : ['TGC','TGT'],
                'Asp' : ['GAT','GAC'],
                'Glu' : ['GAA','GAG'],
                'Phe' : ['TTT','TTC'],
                'Gly' : ['GGC','GGT','GGG','GGA'],
                'His' : ['CAT','CAC'],
                'Ile' : ['ATT','ATC','ATA'],
                'Lys' : ['AAA','AAG'],
                'Leu' : ['CTG','TTG','TTA','CTT','CTC','CTA'],
                'Met' : ['ATG'],
                'Asn' : ['AAC','AAT'],
                'Pro' : ['CCG','CCA','CCT','CCC'],
                'Gln' : ['CAG','CAA'],
                'Arg' : ['CGT','CGC','CGG','CGA','AGA','AGG'],
                'Ser' : ['AGC','TCT','TCC','TCG','AGT','TCA'],
                'Thr' : ['ACC','ACG','ACT','ACA'],
                'Val' : ['GTG','GTT','GTC','GTA'],
                'Trp' : ['TGG'],
                'Tyr' : ['TAT','TAC'],
                ' * ' : ['TAA','TGA','TAG']
                }

        # Search through codon table and find corresponding residue
        translated_raw = []
        for target_codon in codon_split:
            for res, cod in codon_table.items():
                if target_codon in cod:
                    translated_raw.append(res)

        # 
        self.translated_long = ' '.join(translated_raw)            
        return self.translated_long









