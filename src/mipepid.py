import sys

import pandas as pd

import ML
from ORF import ORF


class Dataset:
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

        # read the column labels
        col = infile.readline().rstrip().split(',')

        for line in infile:
            field = line.rstrip().split(',')
            self.data.append({col[i]: field[i] for i in range(len(col))})

        return len(self.data)

    def orf_filter(self, label, constraints):
        """-----------------------------------------------------------------------------------------
        return an list of ORFs object with the ORFs of interest. Constraints can be used to filter
        the list. Constraints is a dictionary. The keys correspond to the column labels and the
        values are arrays that specify allowed values. For "DNAlength" the value indicates the
        minimum acceptable length (including stop codon).

        positive fields: SmProtID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,
        stopCodonSite,transcriptID,corresponding_EnsemblTranscriptIDs,
        corresponding_transcriptBiotypes,transcriptDNAseq,corresponding_EnsemblGeneIDs,
        corresponding_geneBiotypes,dataSource,IsHighConfidence

        negative fields: orfID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,
        EnsemblTranscriptID,transcriptBiotype,transcriptDNAseq,EnsemblGeneID,geneBiotype

        :param self: object         Dataset
        :param label: string        should be "positive" or "negative"
        :param constraints: dict    keys=column labels, values = [] of allowed values
        :return: list of ORF        selected ORFs (see ORF.py)
        -----------------------------------------------------------------------------------------"""
        selected_orfs = []
        if label == "positive":
            for entry in self.data:
                orf = ORF(entry['DNAseq'], seqid=entry['SmProtID'])
                # DNAlength does not include stop codon
                orf.pos = [[0, int(entry['stopCodonSite']) - int(entry['startCodonSite']) + 3]]
                orf.offset = entry['startCodonSite']
                orf.label = 'positive'

                if self.apply_constraints(entry, constraints):
                    selected_orfs.append(orf)

        elif label == "negative":
            for entry in self.data:
                # transcriptDNAseq is the complete sequence, we'll re-call the orfs
                orf = ORF(entry['transcriptDNAseq'], seqid=entry['EnsemblTranscriptID'])
                # DNAlength includes stop codon
                orf.pos = [[0, len(orf.seq) - 3]]
                orf.offset = 0
                orf.label = 'negative'
                orf.split_orfs()

                if self.apply_constraints(entry, constraints):
                    selected_orfs.append(orf)

        else:
            sys.stderr.write(f'mipepid:Dataset:orf_filter() - unkown dataset label ({label})')

        return selected_orfs

    @staticmethod
    def apply_constraints(entry, constraints):
        """-----------------------------------------------------------------------------------------
        Returns true if this orf meets all constraints

        :param entry: dict          information describing a sequence entry
        :param constraints: dict    keys: column titles in entry, values: list of matching terms
        :return: bool               True if constraints are met
        -----------------------------------------------------------------------------------------"""
        ok = True
        for col in constraints:
            if col == 'DNAlength':
                if int(entry['DNAlength']) < constraints[col]:
                    ok = False
                    break
            elif entry[col] not in constraints[col]:
                ok = False
                break

        return ok


# End of class Dataset

# def MiPepid(input_fname, output_fname):
#     """---------------------------------------------------------------------------------------------
#     MiPepid main program
#     Uses Biopython SeqIO to read Fasta sequences from input_fname
#
#     :param input_fname:   string, input file name
#     :param output_fname:  string, output file name
#     :return:
#     ---------------------------------------------------------------------------------------------"""
#     # load the model
#     # logr, threshold = ML.load_model(model_fname='model/newmodel.sav')
#
#     # initialize the output file and write the header
#     columns = ['sORF_ID', 'sORF_seq', 'transcript_DNA_sequence_ID', 'start_at', 'end_at',
#                'classification', 'probability']
#     df = pd.DataFrame(columns=columns)
#     df.to_csv(output_fname, index=False)
#     print('Begin writing the output file: ' + output_fname)
#
#     all_sorfs = []
#     n_predicted = 0
#
    # for rec in SeqIO.parse(input_fname, 'fasta'):
    #     # each ORF sequence is converted to the set of all possible beginning ATGs are the single stop codon. Why?
    #     # Each one of the stop codons could be a coding RF, not necessarily the longes one
    #     orfs = ORFs(str(rec.seq).upper(), rec.id)
    #     all_sORFs += orfs.set()
    #
    #     # Process in batch, each batch >= 1000 ORFs
    #     if len(all_sORFs) > 1000:
    #         ML.batch_predict(all_sORFs, logr, threshold, output_fname)
    #         all_sORFs = []
    #         print(f'Wrote predicted coding/noncoding for {len(all_sORFs)} sORFs')
    #
    # if len(all_sORFs) > 0:
    #     n_predicted += len(all_sORFs)
    #     ML.batch_predict(all_sORFs, logr, threshold, output_fname)
    #     print('Finished writing sORF predictions.')
    #
    # print(f'{n_predicted} sorfs')
    #
    # return n_predicted


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':
    batch_size = 999
    output_fname = 'debug.test.csv'
    model_fname = 'model/newmodel.sav'

    pos = Dataset(sys.argv[1], 'positive')
    filtered_orfs = pos.orf_filter('positive', {})
    neg = Dataset(sys.argv[2], 'negative')
    filtered_orfs += neg.orf_filter('negative', {'DNAlength': 50})

    # load the model
    logr, threshold = ML.load_model(model_fname=model_fname)

    batch = []
    n_predicted = 0
    for s in filtered_orfs:

        if len(s.pos) > 1:
            # multiple ORFs
            n_orf = 0
            for begin, end in s.pos:
                n_orf += 1
                sorf_seq = s.seq[begin:end]
                sorf_id = s.seqid + f'_ORF{n_orf:02d}'
                batch.append([sorf_id, sorf_seq, sorf_id, begin, end])
        else:
            # just one ORF
            batch.append([s.seqid, s.seq, s.seqid, s.pos[0][0], s.pos[0][1]])

        # Process in batch, each batch >= batch_size ORFs
        if len(batch) > batch_size:
            ML.batch_predict(batch, logr, threshold, output_fname)
            print(f'\t... Predicting coding/noncoding for {len(batch)} sORFs/ total={n_predicted}')
            n_predicted += len(batch)
            batch = []

    if len(batch) > 0:
        n_predicted += len(batch)
        ML.batch_predict(batch, logr, threshold, output_fname)
        n_predicted += len(batch)
        print(f'\t... Predicting coding/noncoding for {len(batch)} sORFs/ total={n_predicted}')

    print(f'{n_predicted} sORF predictions written to {output_fname}')

exit(0)
