import sys

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
                orf.tag = Dataset.to_unique_tag(entry['corresponding_transcriptBiotypes'])
                if entry['IsHighConfidence'] == 'Yes':
                    orf.tag.append('HighConfidence')

                if self.apply_constraints(entry, constraints):
                    selected_orfs.append(orf)

        elif label == "negative":
            negative_id_list = {}
            for entry in self.data:
                if entry['EnsemblTranscriptID'] in negative_id_list:
                    # skip any duplicate sequence entries
                    continue
                negative_id_list[entry['EnsemblTranscriptID']] = 1

                # transcriptDNAseq is the complete sequence, we'll re-call the orfs
                orf = ORF(entry['transcriptDNAseq'], seqid=entry['EnsemblTranscriptID'])
                # DNAlength does not include stop codon
                orf.pos = [[0, len(orf.seq) - 3]]
                orf.offset = 0
                orf.label = 'negative'
                orf.tag = Dataset.to_unique_tag(entry['transcriptBiotype'])
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

    @staticmethod
    def to_unique_tag(tag, separator=';'):
        """-----------------------------------------------------------------------------------------
        Convert s string of multiple tags separated by semicolons to a unique list
        protein_coding;protein_coding;junk -> ['protein_coding', 'junk']

        :param tag: string      one or more tags, delimited by separator
        :param separator:       delimiter for fields in string
        :return: list           unique list of tag strings
        -----------------------------------------------------------------------------------------"""
        field = tag.strip().split(separator)
        tagdict = {}
        for t in field:
            try:
                tagdict[t] += 1
            except KeyError:
                tagdict[t] = 1

        taglist = list(tagdict)
        if not taglist:
            taglist = []

        return taglist


# ===================================================================================================
# End of class Dataset
# ===================================================================================================


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':
    batch_size = 999
    output_fname = 'debug.test.csv'
    model_fname = 'model/newmodel.sav'
    k = 4
    filter_len = 30
    threshold = 0.75

    pos = Dataset(sys.argv[1], 'positive')
    sys.stderr.write(f'Positive sequences read: {len(pos.data)}\n')
    filtered_orfs = pos.orf_filter('positive', {'DNAlength': filter_len})
    n_pos_filt = len(filtered_orfs)
    sys.stderr.write(f'Positive sequences after filtering ({filter_len}): {n_pos_filt}\n')

    neg = Dataset(sys.argv[2], 'negative')
    sys.stderr.write(f'Negative sequences read: {len(neg.data)}\n')
    filtered_orfs += neg.orf_filter('negative', {'DNAlength': filter_len})
    n_neg_filt = len(filtered_orfs) - n_pos_filt
    sys.stderr.write(f'Negative sequences after filtering ({filter_len}): {n_neg_filt}\n')

    # load the model
    sys.stderr.write(f'\nLoading regression model from {model_fname}...\n')
    logr, threshold = ML.load_model(model_fname=model_fname)
    sys.stderr.write(f'... completed loading model\n')

    batch = []
    n_predicted = 0
    sys.stderr.write(f'\nBeginning predictions for {len(filtered_orfs)} sequences\n')
    for s in filtered_orfs:

        sorf_tag = ';'.join(s.tag)
        # if s.seqid.startswith('ENST'):
        #     print('test')
        if len(s.pos) > 1:
            # multiple ORFs
            n_orf = 0
            for begin, end in s.pos:
                n_orf += 1
                sorf_seq = s.seq[begin:end]
                if len(sorf_seq) < k:
                    # if sequence is less than k, logistic regression fails (no features)
                    continue
                sorf_id = s.seqid + f'_ORF{n_orf:02d}'
                batch.append([sorf_id, sorf_seq, begin, end, s.label, sorf_tag])
        else:
            # just one ORF
            begin = s.pos[0][0]
            end = s.pos[0][1]
            batch.append([s.seqid, s.seq[begin:end], begin, end, s.label, sorf_tag])

        # Process in batch, each batch >= batch_size ORFs
        if len(batch) > batch_size:
            ML.batch_predict(batch, logr, threshold, output_fname)
            n_predicted += len(batch)
            print(f'... Predicting coding/noncoding for {len(batch)} sORFs: total={n_predicted}')
            batch = []

    if len(batch) > 0:
        n_predicted += len(batch)
        ML.batch_predict(batch, logr, threshold, output_fname)
        n_predicted += len(batch)
        print(f'... Predicting coding/noncoding for {len(batch)} sORFs: total={n_predicted}')

    print(f'\n{n_predicted} sORF predictions written to {output_fname}')

exit(0)
