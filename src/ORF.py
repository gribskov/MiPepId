class ORF:
    """=============================================================================================
    An orf holds a sequence, a list of begin and end coordinates of ORFs, and minimal metadata
    sequence ID, source ID

    sequence    entire subject sequence (preferably) or an extracted subsequence defined by offset
    offset      offset that would convert internal coordinate to the source coordinate
    seqid       ID of subject sequence
    label       text label describing the sequence ('positive', 'negative', 'unknown' ...)
    tag         additional class tags such as 'high quality' or 'miRNA'
    pos: list   [begin,end] positions of ORFs beginning with start and ending with stop codons

    ============================================================================================="""
    allowed_start_codons = ['ATG']
    allowed_stop_codons = ['TAA', 'TAG', 'TGA']

    def __init__(self, seq, seqid='', start=None, stop=None, split=False):
        """-----------------------------------------------------------------------------------------
        attributes:
        seq: string             sequence of the region of interest, preferably entire source sequence
        id: string              ID of source sequence
        pos: list of list       [[orf_begin_pos, orf_end_pos], ... ]
        offset: int             offset that would convert internal coordinate to the source 
                                coordinate
        label: string           label, expect 'positive' or 'negative'
        tag: list of string     additional class tags such as 'high quality' or 'miRNA'

        :param seq: string      populates ORF.seq
        :param seqid: string    populates ORF.seqid
        :param split: bool      if True, the provided region will be split into smaller ORFs
        -----------------------------------------------------------------------------------------"""
        self.seq = seq
        self.seqid = ''
        self.pos = [[0, len(self.seq)]]
        self.offset = 0
        self.label = ''
        self.tag = []

        # override defaults
        if id:
            self.seqid = seqid
        if start:
            self.start = start
        if stop:
            self.stop = stop
        if split:
            self.split_orfs()

    def split_orfs(self):
        """-----------------------------------------------------------------------------------------
        search the current ORF regions stored in pos and find minimal ORFs that start at a start
        codon and end at the first encountered stop codon.

        Updates pos attribute in object. Positions give the position of the first base of the
        start/stop codon in the sequence

        TODO should terminal ORFs without stop codon be reported?

        :return: int    number of ORFs found
        -----------------------------------------------------------------------------------------"""
        sequence = self.seq
        start = ORF.allowed_start_codons
        stop = ORF.allowed_stop_codons
        rflist = []

        # orfopen gives the current start position in each RF. ORF must begin with a start codon
        orfopen = [-1, -1, -1]
        for pos in range(len(sequence)):
            frame = pos % 3
            codon = sequence[pos:pos + 3]
            if codon in start:
                # start codon
                if orfopen[frame] == -1:
                    # no current working orf in this frame, start a new one
                    orfopen[frame] = pos

            elif codon in stop:
                # stop codon
                if orfopen[frame] > -1:
                    # there is a current orf in this frame, save and close
                    rflist.append([orfopen[frame], pos])
                    orfopen[frame] = -1

        self.pos = rflist
        return len(rflist)

    def long_orfs(self):
        """-----------------------------------------------------------------------------------------
        for each reading frame, remove all in-frame start and stop codons creating a long ORF

        :return: list   three sequence strings
        -----------------------------------------------------------------------------------------"""
        seq = self.seq
        orf = []
        for frame in (0, 1, 2):
            newseq = ''
            for pos in range(frame, len(seq)-3, 3):
                codon = seq[pos:pos+3]
                if codon in self.start or codon in self.stop:
                    continue
                newseq += codon
            orf.append(newseq)

        return orf

# ==================================================================================================
# end of class ORF
# ==================================================================================================

# testing

if __name__ == '__main__':
    seq = ORF('AACTGCCCTCAAGGAGCTTACAATCTAGCTGGGGGTAAATGACTTGCACATGAACACAACTAGACTGTGAGCTTCTAGAGGGC')
    orf = seq.long_orfs()

    exit(0)
