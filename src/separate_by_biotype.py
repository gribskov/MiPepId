"""=================================================================================================
separate the original negative csv file into fasta files by biotype. fields:
orfID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,EnsemblTranscriptID,transcriptBiotype,transcriptDNAseq,EnsemblGeneID,geneBiotype

we want orfID (field 0), DNAseq (field 1), transcriptBiotype (field 8)

Michael Gribskov     25 August 2023
================================================================================================="""
import sys

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    csv = open(sys.argv[1], 'r')

    btype = {}
    header = csv.readline();
    seq_n = 0
    for line in csv:
        seq_n += 1
        field = line.strip().split(',')
        id = field[0]
        seq = field[1]
        biotype = field[8]

        # print( f'{seq_n}\t{id}\t{biotype}\n{seq}')

        if biotype in btype:
            btype[biotype].append( {'id':id, 'seq':seq, 'biotype':biotype})
        else:
            btype[biotype] = [{'id':id, 'seq':seq, 'biotype':biotype}]

    print(f'\n{seq_n} ORFs read')
    print(f'\nSeparating by biotype')
    for biotype in btype:
        print(f'{biotype}\t{len(btype[biotype])}')
        out = open(biotype+'.fa', 'w')
        for entry in btype[biotype]:
            out.write(f">{entry['id']}\t{biotype}\n{entry['seq']}\n")

        out.close()

    exit(0)
