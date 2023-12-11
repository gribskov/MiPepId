class ORF:
    """=================================================================================================================
    ORF holds the sequence of a single ORF, and its position in the original sequence
    ================================================================================================================="""

    def __init__(self, seq, start_codon_site=None, stop_codon_site=None):
        self.seq = seq
        self.start_codon = seq[:3]
        self.stop_codon = seq[-3:]
        self.length = len(seq)
        self.start_codon_site = start_codon_site
        self.stop_codon_site = stop_codon_site


class ORFs:
    """=================================================================================================================

    ================================================================================================================="""

    def __init__(self, DNA_sequence, candidate_start_codons=['ATG'], candidate_stop_codons=['TAA', 'TAG', 'TGA']):

        self.seq = DNA_sequence
        self.candidate_start_codons = candidate_start_codons
        self.candidate_stop_codons = candidate_stop_codons

        self.all_ORFs = []

        for i in range(3):  # the 3 frames
            S = self.seq[i:]
            fragment_regions = self.__break_sequence_into_fragments_by_stopCodon__(S)

            if len(fragment_regions) > 0:
                for fragment_region in fragment_regions:
                    start_codon_starting_sites = self.__find_the_starting_sites_of_all_startCodons__(fragment_region, S)
                    if len(start_codon_starting_sites) > 0:
                        for elm in start_codon_starting_sites:
                            start_codon_site = i + elm;
                            stop_codon_site = i + fragment_region[1] - 3
                            seq = DNA_sequence[start_codon_site:(stop_codon_site + 3)]
                            obj = ORF(seq, start_codon_site, stop_codon_site)
                            self.all_ORFs.append(obj)

    def get_orfs_original(self):
        """-----------------------------------------------------------------------------------------
        Mengmengs original code for getting the Start-Stop ORFs in the sequence. Originally placed
        in __init__()
        -----------------------------------------------------------------------------------------"""
        for i in range(3):  # the 3 frames
            fragment_regions = self.__break_sequence_into_fragments_by_stopCodon__(i)

            if len(fragment_regions) > 0:
                for fragment_region in fragment_regions:
                    start_codon_starting_sites = self.__find_the_starting_sites_of_all_startCodons__(fragment_region, S)
                    if len(start_codon_starting_sites) > 0:
                        for elm in start_codon_starting_sites:
                            start_codon_site = i + elm;
                            stop_codon_site = i + fragment_region[1] - 3
                            seq = DNA_sequence[start_codon_site:(stop_codon_site + 3)]
                            obj = ORF(seq, start_codon_site, stop_codon_site)
                            self.all_ORFs.append(obj)

    def get_orfs(self):
        """---------------------------------------------------------------------------------------------
        return a non-overlapping list of orfs beginning with any codon in start and ending with one of
        the infram stop codons in stop

        :param sequence:    string with complete sequence
        :param start:       list of start codons, e.g. ['ATG']
        :param stop:        list of stop codons, e.g. ['TAA', 'TAG', 'TGA']
        :return:            list of [ORF_string, begin, end]
        ---------------------------------------------------------------------------------------------"""
        sequence = self.DNA_sequence
        start = self.candidate_start_codons
        stop = self.candidate_stop_codons
        rflist = []
        open = [False, False, False]
        for pos in range(len(sequence)):
            frame = pos % 3
            if sequence[pos:pos + 3] in start:
                # start codon
                if not open[frame]:
                    # no current working orf in this frame, start a new one
                    open[frame] = pos

            elif sequence[pos:pos + 3] in stop:
                # stop codon
                if open[frame]:
                    # there is a current orf in this frame, save and close
                    rflist.append([sequence[open[frame]:pos + 3], open[frame], pos + 2])
                    open[frame] = False

        return rflist

    def __break_sequence_into_fragments_by_stopCodon__(self, frame):
        """-------------------------------------------------------------------------------------------------------------
        Given a DNA sequence, break this into fragments with each fragment ending with a stop codon and the length of
        the fragment is a multiple of 3. Only the first frame of this DNA sequence is considered.
       mrg: reworked to not copy sequence
        Returns:
          a list of the regions of the fragments. Each region in the list is a tuple consisting of 2 elements: the
          starting site of the fragment (0-based) and the end site of the fragment (1-based), so that the 2 numbers can
          be used directly in Python slicing.
        -------------------------------------------------------------------------------------------------------------"""
        stop_codon_starting_sites = []
        for j in range(frame, len(self.seq)-3, 3):
            if self.seq[j:j + 3] in self.candidate_stop_codons:
                stop_codon_starting_sites.append(j)

        fragment_regions = []

        if len(stop_codon_starting_sites) > 0:
            start = 0
            for elt in stop_codon_starting_sites:
                this_region = (start, elt + 3)
                fragment_regions.append(this_region)
                start = elt + 3

        return fragment_regions

    def __find_the_starting_sites_of_all_startCodons__(self, fragment_region, S):
        """-------------------------------------------------------------------------------------------------------------
        Given a DNA sequence and the region for a fragment, find the starting site of every start codon in the first
        translating frame of this fragment.
        -------------------------------------------------------------------------------------------------------------"""
        start_codon_starting_sites = []
        for k in range(fragment_region[0], fragment_region[1] - 3, 3):
            if S[k:(k + 3)] in self.candidate_start_codons:
                start_codon_starting_sites.append(k)
        return start_codon_starting_sites


# End of class ORFs


def collect_and_name_sORFs_from_an_ORFs_object(obj_ORFs, transcript_seq_ID):
    """-----------------------------------------------------------------------------------------------------------------
    Creates a list of all possible ATG -> stop ORFs from a list of maximal ORFS. These sORFs overlap
    :param obj_ORFs:
    :param transcript_seq_ID:
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    sORFs = []
    count = 0
    for ORF in obj_ORFs.all_ORFs:
        if ORF.length <= 303:
            count += 1
            orfID = transcript_seq_ID + '_ORF' + str(count)
            this_sORF = [orfID, ORF.seq, transcript_seq_ID, ORF.start_codon_site + 1, ORF.stop_codon_site + 3]
            sORFs.append(this_sORF)

    return sORFs


def get_orfs_simple(sequence, start=['ATG'], stop=['TAA', 'TAG', 'TGA']):
    """---------------------------------------------------------------------------------------------
    return a non-overlapping list of orfs beginning with any codon in start and ending with one of
    the infram stop codons in stop

    :param sequence:    string with complete sequence
    :param start:       list of start codons, e.g. ['ATG']
    :param stop:        list of stop codons, e.g. ['TAA', 'TAG', 'TGA']
    :return:            list of [ORF_string, begin, end]
    ---------------------------------------------------------------------------------------------"""
    rflist = []
    open = [False, False, False]
    for pos in range(len(sequence)):
        frame = pos % 3
        if sequence[pos:pos + 3] in start:
            # start codon
            if not open[frame]:
                # no current working orf in this frame, start a new one
                open[frame] = pos

        elif sequence[pos:pos + 3] in stop:
            # stop codon
            if open[frame]:
                # there is a current orf in this frame, save and close
                rflist.append([sequence[open[frame]:pos + 3], open[frame], pos])
                open[frame] = False

    return rflist
