"""=================================================================================================
data in the negative_original_data.csv and file contains some overlapping
orfs. Regenerate the data from the source transcript DNA sequences in the file.

Michael Gribskov     07 December 2023
================================================================================================="""
import sys
import pandas
from ORF import get_orfs_simple


def write_data(file, df):
    """---------------------------------------------------------------------------------------------
    write the list of ORFs as CSV. Columns are
    orfID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,EnsemblTranscriptID,
    transcriptBiotype,transcriptDNAseq,EnsemblGeneID,geneBiotype

    :param file: filepointer open for writing
    :param df:   dataframe from original csv with ORF column added
    :return:     rows written
    ---------------------------------------------------------------------------------------------"""
    nseq = 0
    nrow = 0
    file.write('orfID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,EnsemblTranscriptID,transcriptBiotype,transcriptDNAseq,EnsemblGeneID,geneBiotype\n')
    for index, row in df.iterrows():
        nseq += 1
        norf = 0
        for orf in row['ORF']:
            nrow += 1
            norf += 1
            orfid = f"{row['EnsemblTranscriptID']}_ORF{norf:02d}"
            file.write('{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(
                orfid,
                orf[0],
                row['DNAlength'],
                orf[1],
                orf[2],
                orf[0][:3],
                orf[0][-3:],
                row['EnsemblTranscriptID'],
                row['transcriptBiotype'],
                row['transcriptDNAseq'],
                row['EnsemblGeneID'],
                row['geneBiotype']
                ))

    return nrow

def write_fa(file, df):
    """---------------------------------------------------------------------------------------------
    Write the orf out in fasta format

    :param file: filehandle open for reading
    :param df:   dataframe with sorf data
    :return:     number of sequences written
    ---------------------------------------------------------------------------------------------"""
    nseq = 0
    for index, row in df.iterrows():
        nseq += 1
        file.write(f'>{row["orfID"]}\n{row["DNAseq"]}\n')

    return nseq

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    orig = pandas.read_csv(sys.argv[1])
    print(f'{len(orig)} rows read from {sys.argv[1]}')

    # orf01 = orig[orig['orfID'].str.endswith('ORF01')]
    orf01 = pandas.DataFrame(orig.loc[orig['orfID'].str.endswith('ORF01'),])
    print(f'{len(orf01)} sequences found')
    orf01['ORF'] = orf01['transcriptDNAseq'].apply(get_orfs_simple)

    out = open(sys.argv[2], 'w')
    nrow = write_data(out, orf01)
    out.close()
    print(f'{nrow} rows written to {sys.argv[2]}')

    faname = sys.argv[2] + '.fa'
    fa_out = open(faname, 'w')
    nrows = write_fa(fa_out, orf01)
    print(f'{nrow} rows written to {faname}')

    exit(0)
