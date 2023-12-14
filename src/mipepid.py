from ORF import ORF, ORFs, collect_and_name_sORFs_from_an_ORFs_object
import ML

from Bio import SeqIO
import pandas as pd
import sys


class Dataset():
    """=============================================================================================
    Mengmeng's positive and negative datasets have different fields, so the simplest solution is to
    store them as a dict and use functions to convert them to ORFs for prediction.

    the positive examples include only the ORF of interest in the DNASeq field:
        fields: SmProtID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,
        transcriptID,corresponding_EnsemblTranscriptIDs,corresponding_transcriptBiotypes,
        transcriptDNAseq,corresponding_EnsemblGeneIDs,corresponding_geneBiotypes,dataSource,
        IsHighConfidence

    the negative examples contain multiple ORFs from the same sequence: while reading only the first
    transcript ID is stored. The ORFs are called later when adding to ORF lists.
        fields: orfID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,
        EnsemblTranscriptID,transcriptBiotype,transcriptDNAseq,EnsemblGeneID,geneBiotype

    a Dataset is a list of dictionaries with the above fields
    ============================================================================================="""

    def __init__(self, filename, label):
        """-----------------------------------------------------------------------------------------
        Constructor

        :param label: string        expect "positive" or "negative"
        :param filename: string     filename for the input data
        -----------------------------------------------------------------------------------------"""
        self.data = []
        if label:
            self.label = label

        if filename:
            self.read_data(filename)

    def read_data(self, filename):
        """-----------------------------------------------------------------------------------------
        Read in the data and store as annotated in the .csv file as a list of hashes. The
        first line is a comma delimited string with the column titles. These become the  dictionary
        keys.

        :param filename: string     path to the source .csv file
        :return: int                number of sequences read
        -----------------------------------------------------------------------------------------"""
        infile = open(filename, 'r')
        if not infile:
            sys.stderr.write(f'mipepid.py:read_data() cannot open input file ({filename}')
            exit(1)

        col = infile.readline().rstrip().split(',')

        n_read = 0
        for line in infile:
            field = line.rstrip().split(',')
            self.data.append({col[i]: field[i] for i in range(len(col))})

        return len(self.data)

    def orf_filter(self, label, constraints):
        """-----------------------------------------------------------------------------------------
        return an ORFs object with a list of the ORFs of interest. Constraints can be used to filter
        the list. Constraints is a dictionary. The keys correspond to the column labels and the
        values are arrays that specify allowed values. For "DNAlength" the value indicates the
        minimum acceptable length (including stop codon).

        SmProtID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,
        transcriptID,corresponding_EnsemblTranscriptIDs,corresponding_transcriptBiotypes,
        transcriptDNAseq,corresponding_EnsemblGeneIDs,corresponding_geneBiotypes,dataSource,
        IsHighConfidence

        :param self: object         Dataset
        :param label: string        should be "positive" or "negative"
        :param constraints: dict    keys=column labels, values = [] of allowed values
        :return: object             ORFs (see ORF.py)
        -----------------------------------------------------------------------------------------"""
        selected_orfs  = ORFs()
        if label == "positive":
            for entry in self.data:
                seq = entry['DNAseq']
                orf = ORF(seq, entry['startCodonSite'], entry['stopCodonSite'])
                orf.start_codon = seq[:3]
                orf.stop_codon = seq[-3:]
                if apply_constraints( orf, entry, constraints):
                    selected_orfs.orflist.append(orf)

        elif label == "negative":
            for entry in self.data:

        else:
            sys.stderr.write(f'mipepid:Dataset:orf_filter() - unkown dataset label ({label})')

        return len(orflist)

    @staticmethod
    def allply_constraints(entry, constraints):
        """-----------------------------------------------------------------------------------------
        Returns true if this orf meets all constraints

        :param entry:
        :param constraints:
        :return:
        -----------------------------------------------------------------------------------------"""
        OK = True

        return OK


# End of class Dataset

def MiPepid(input_fname, output_fname):
    """---------------------------------------------------------------------------------------------
    MiPepid main program
    Uses Biopython SeqIO to read Fasta sequences from input_fname

    :param input_fname:   string, input file name
    :param output_fname:  string, output file name
    :return:
    ---------------------------------------------------------------------------------------------"""
    # initialize the output file and write the header
    columns = ['sORF_ID', 'sORF_seq', 'transcript_DNA_sequence_ID', 'start_at', 'end_at',
               'classification', 'probability']
    df = pd.DataFrame(columns=columns)
    df.to_csv(output_fname, index=False)
    print('Begin writing the output file: ' + output_fname)

    # load the model
    logr, threshold = ML.load_model(model_fname='model/newmodel.sav')

    # gather all the sORFs and write
    all_sORFs = []
    n_predicted = 0

    for rec in SeqIO.parse(input_fname, 'fasta'):
        # each ORF sequence is converted to the set of all possible beginning ATGs are the single stop codon. Why?
        # Each one of the stop codons could be a coding RF, not necessarily the longes one
        orfs = ORFs(str(rec.seq).upper(), rec.id)
        all_sORFs += orfs.set()

        # Process in batch, each batch >= 1000 ORFs
        if len(all_sORFs) > 1000:
            ML.predict_on_one_batch_and_write(all_sORFs, logr, threshold, output_fname)
            all_sORFs = []
            print(f'Wrote predicted coding/noncoding for {len(all_sORFs)} sORFs')

    if len(all_sORFs) > 0:
        n_predicted += len(all_sORFs)
        ML.predict_on_one_batch_and_write(all_sORFs, logr, threshold, output_fname)
        print('Finished writing sORF predictions.')

    print(f'{n_predicted} sorfs')

    return n_predicted


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':
    input_fname = output_fname = None

    pos = Dataset(sys.argv[1], 'positive')

    if len(sys.argv) == 1:
        print('Please specifiy the input fasta file using the following command:')
        print('python3 ./src/mipepid.py input_file_name')
    elif len(sys.argv) == 2:
        input_fname = sys.argv[1]
        output_fname = 'Mipepid_results.csv'
    else:
        input_fname = sys.argv[1]
        output_fname = sys.argv[2]

    if input_fname and output_fname:
        MiPepid(input_fname, output_fname)
