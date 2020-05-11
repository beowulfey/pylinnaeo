from abc import ABC

from Bio import SeqRecord


class SeqR(SeqRecord.SeqRecord, ABC):
    """
    Implementation of SeqRecord that allows for comparisons.
    ID is fixed. It is the starting name. If a fasta is properly formatted, it
    should be the ID number of the protein.
    Otherwise, the ID and name end up the same.
    However, the name is the part that is always changed, and does not come with exporting.
    # TODO: exporting sequence, allow for choosing fasta format, or name only!
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

    #def sName(self):
    #    return self._sname

    #def setSeqName(self, name):
    #    self._sname = name

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
        if self.seq == other.seq:
            return True
        else:
            return False
