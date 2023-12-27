"""=================================================================================================
MiPepid:long_negative.py

generate long negative sequences by removing in-frame start and stop codons in each RF

27 December 2023     gribskov
================================================================================================="""
import sys
import pandas
from ORF import ORF

# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    neg = pandas.read_csv(sys.argv[1])
    rf = [[], [], []]
    ldist = []
    for index, entry in neg.iterrows():
        if not entry['orfID'].endswith('ORF01'):
            continue

        s = ORF(seq=entry['transcriptDNAseq'])
        frame = s.long_orfs()
        for f in range(0, 3):
            this_orf = entry.copy(deep=True)
            this_orf['orfID'] = this_orf['orfID'].replace('_ORF01', f'_LRF{f}')
            this_orf['DNAseq'] = frame[f]
            this_orf['DNAlength'] = len(frame[f])
            this_orf['startCodon'] = frame[f][0:3]
            this_orf['stopCodon'] = frame[f][-3:]
            this_orf['startCodonSite'] = 0
            this_orf['stopCodonSite'] = len(frame[f])  # stop codon is after end of sequence
            rf[f].append(this_orf)
            ldist.append({'id':this_orf['orfID'], 'len':len(frame[f])})

    for f in range(0, 3):
        df = pandas.DataFrame(rf[f])
        lrf = f'negative_lrf_{f}'
        df.to_csv(lrf)

    # write out lengths
    lout = open('length_out.txt', 'w')
    lout.write('id\tlength\n')
    for row in ldist:
        lout.write(f"{row['id']}\t{row['len']}\n")
    lout.close()

    exit(0)
