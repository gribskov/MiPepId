"""=================================================================================================


Michael Gribskov     07 December 2023
================================================================================================="""


def get_orfs(sequence, start, stop):
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
                rflist.append([sequence[open[frame]:pos + 3], open[frame], pos + 2])
                open[frame] = False

    return rflist


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    seqs = ['CGCCGGCCGATGGGCGTCTTACCAGACATGGTTAGACCTGGCCCTCTGTCTAATACTGTCTGGTAAAACCGTCCATCCGCTGC',
            'CTGTGTGTGATGAGCTGGCAGTGTATTGTTAGCTGGTTGAATATGTGAATGGCATCGGCTAACATGCAACTGCTGTCTTATTGCATATACA',
            'ACCCAAACCCTAGGTCTGCTGACTCCTAGTCCAGGGCTCGTGATGGCTGGTGGGCCCTGAACGAGGGGTCTGGAGGCCTGGGTTTGAATATCGACAGC',
            'CTGGGGTACGGGGATGGATGGTCGACCAGTTGGAAAGTAATTGTTTCTAATGTACTTCACCTGGTCCACTAGCCGTCCGTATCCGCTGCAG',
            'TGGGAAACATACTTCTTTATATGCCCATATGGACCTGCTAAGCTATGGAATGTAAAGAAGTATGTATCTCA',
            'TACTTAAAGCGAGGTTGCCCTTTGTATATTCGGTTTATTGACATGGAATATACAAGGGCAAGCTCTCTGTGAGTA',
            'CCGCCCCGGGCCGCGGCTCCTGATTGTCCAAACGCAATTCTCGAGTCTATGGCTCCGGCCGAGAGTTGAGTCTGGACGTCCCGAGCCGCCGCCCCCAAACCTCGAGCGGG',
            'CTGACTATGCCTCCCCGCATCCCCTAGGGCATTGGTGTAAAGCTGGAGACCCACTGCCCCAGGTGCTGCTGGGGGTTGTAGTC',
            'CACTCTGCTGTGGCCTATGGCTTTTCATTCCTATGTGATTGCTGTCCCAAACTCATGTAGGGCTAAAAGCCATGGGCTACAGTGAGGGGCGAGCTCC',
            'ACTGTCCTTTTTCGGTTATCATGGTACCGATGCTGTATATCTGAAAGGTACAGTACTGTGATAACTGAAGAATGGTGGT',
            'CCTCAGAAGAAAGATGCCCCCTGCTCTGGCTGGTCAAACGGAACCAAGTCCGTCTTCCTGAGAGGTTTGGTCCCCTTCAACCAGCTACAGCAGGGCTGGCAATGCCCAGTCCTTGGAGA',
            'AGAGATGGTAGACTATGGAACGTAGGCGTTATGATTTCTGACCTATGTAACATGGTCCACTAACTCT',
            'GGTCTCTGTGTTGGGCGTCTGTCTGCCCGCATGCCTGCCTCTCTGTTGCTCTGAAGGAGGCAGGGGCTGGGCCTGCAGCTGCCTGGGCAGAGCGG',
            'AAACGATACTAAACTGTTTTTGCGATGTGTTCCTAATATGCACTATAAATATATTGGGAACATTTTGCATGTATAGTTTTGTATCAATATA']

    start = ['ATG']
    stop = ['TAA', 'TAG', 'TGA']

    for s in seqs:
        print(s)
        orfs = get_orfs(s, start, stop)
        for rf in orfs:
            print(f'    {rf[1]:4d}    {rf[2]:4d}    {rf[0]}')

        print()

    exit(0)
