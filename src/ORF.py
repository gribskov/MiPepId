class ORF:
    """=============================================================================================
    An orf holds a sequence, a list of begin and end coordinates of ORFs, and minimal metadata
    sequence ID, source ID
    ============================================================================================="""
    allowed_start_codons = ['ATG']
    allowed_stop_codons = ['TAA', 'TAG', 'TGA']

    def __init__(self, seq, seqid='', start=None, stop=None, split=False):
        """-----------------------------------------------------------------------------------------
        seq: string             sequence of the region of interest, preferably the entire source 
                                sequence
        id: string              ID of source sequence
        start: list of string   allowed start codons, default is 
        stop: list of string    allowed stop codons, default is 
        pos: list of list       [[orf_begin_pos, orf_end_pos], ... ]
        offset: int             offset that would convert internal coordinate to the source 
                                coordinate
        #TODO need to add source sequence ID
        split: bool             if True, the provided region will be split into smaller ORFs 
        -----------------------------------------------------------------------------------------"""
        self.seq = seq
        self.seqid = ''
        self.start = ORF.allowed_start_codons
        self.stop = ORF.allowed_stop_codons
        self.pos = [[0, len(self.seq)]]
        self.offset = 0
        self.label = ''

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

        :return: bool   True
        -----------------------------------------------------------------------------------------"""
        sequence = self.seq
        start = self.start
        stop = self.stop
        rflist = []
        # open gives the current start position in each RF. Must begin with a start codon
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
        return True


# ===================================================================================================
# end of class ORF
# ===================================================================================================
