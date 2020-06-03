from abc import ABC

from Bio import SeqRecord


class SeqR(SeqRecord.SeqRecord, ABC):
    """
    Implementation of SeqRecord that allows for comparisons.
    ID is fixed. It is the starting name. If a fasta is properly formatted, it
    should be the ID number of the protein.
    Otherwise, the ID and name end up the same.
    However, the name is the part that is always changed, and does not show up when exporting.
    # TODO: exporting sequence, allow for choosing fasta format, or name only!
    # TODO: Fix "format" option so that it doesn't include "unknown description"
    """
    def __init__(self,
                 seq,
                 id="<unknown id>",
                 name="<unknown name>",
                 description="<unknown description>",
                 dbxrefs=None,
                 features=None,
                 annotations=None,
                 letter_annotations=None,
                 ):
        super(SeqR, self).__init__(seq, id, name, description,
                                   dbxrefs, features, annotations, letter_annotations)

    def convert(self, seqr):
        """ Converts a Bio.Seqr into my own SeqR. Skips the seq (should be already there)"""
        self.id = seqr.id
        self.name=seqr.id
        self.description = seqr.description
        self.dbxrefs=seqr.dbxrefs
        self.features=seqr.features
        self.annotations=seqr.annotations
        self.letter_annotations=seqr.letter_annotations

    def __lt__(self, other):
        """Define the less-than operand."""
        if self.seq < other.seq:
            return True
        else:
            return False

    def __gt__(self, other):
        """Define the greater-than operand."""
        if self.seq > other.seq:
            return True
        else:
            return False

    def __eq__(self, other):
        """Define the equal-to operand."""
        if type(other) == SeqR:
            if self.seq == other.seq:
                return True
            else:
                return False
        elif type(other) == list:
            for n in other:
                if self.seq == n.seq:
                    return True
                else:
                    return False
