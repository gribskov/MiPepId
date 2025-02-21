"""


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

    unknowns are read in FastA format by and separated into ORFs read_fasta()
"""
# TODO update and incorporate this into current documentation

import sys

import ML
from ORF import ORF


class DataEntry:
    """=============================================================================================
    DataEntry is ultimately used to feed ML:batch_predict: All Data entries have the following
    necessary information for each ORF in the dataset:
        orf_id      name after extracting from original sequence
        orf_seq     sequence after extracting from sequence
        source:     id of original sequence the ORF is extracted from
        start_at    position of the first base of the start codon in the original sequence
        end_at      position of the first base of the stop codon in the original sequence
        class_label typically 'positive', 'negative', or 'unknown'
        tags        list of additional tags (such as data source) to display in report

    optional additional information is stored in
        data        other information from input stored as dict

    Only ORF_seq is used for the prediction; the other fields are just for the report
    ============================================================================================="""

    def __init__(self, orf_id='', orf_seq='', source='', start_at=None, end_at=None,
                 class_label='unknown'):
        """-----------------------------------------------------------------------------------------
        Constructor: see class DataEntry docstring
        TODO change data to a more informative name. info? attributes (like GFF)?
        -----------------------------------------------------------------------------------------"""
        self.orf_id = orf_id
        self.orf_seq = orf_seq
        self.source = source
        self.start_at = start_at
        self.end_at = end_at
        self.class_label = class_label
        self.tags = []
        self.data = []


class DataSet(list):
    """=============================================================================================
    DataSet is a collection of DataEntry
    ============================================================================================="""

    def __init__(self, filename, label):
        """-----------------------------------------------------------------------------------------
        Constructor: see class docstring
        TODO change data to a more informative name. info? attributes (like GFF)?
        -----------------------------------------------------------------------------------------"""
        super().__init__()      # doesn't do anything
        if label:
            self.label = label

        if filename:
            if self.label == 'positive':
                self.read_csv(filename)
            elif self.label == 'negative':
                self.read_csv(filename)
            else:
                # anything else is unknown
                self.read_fasta(filename)

    def read_csv(self, filename):
        """-----------------------------------------------------------------------------------------
        Positive and negative data are specially formatted with metadata that comes from the
        respective sources, SmProt and Ensembl. The fields for each are shown below

        positive
        SmProtID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,transcriptID,
        corresponding_EnsemblTranscriptIDs,corresponding_transcriptBiotypes,transcriptDNAseq,
        corresponding_EnsemblGeneIDs,corresponding_geneBiotypes,dataSource,IsHighConfidence

        negative
        orfID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,
        EnsemblTranscriptID,transcriptBiotype,transcriptDNAseq,EnsemblGeneID,geneBiotype

        :param filename: string     path to the source .csv file
        :return: int                number of sequences read
        -----------------------------------------------------------------------------------------"""
        infile = open(filename, 'r')
        if not infile:
            sys.stderr.write(f'DataSet:read_csv() cannot open input file ({filename}')
            exit(1)

        # read the column labels
        col = infile.readline().rstrip().split(',')
        colidx = {col[i]: i for i in range(len(col))}

        # select columns to populate fields in DataEntry
        if self.label == 'positive':
            seq = colidx['DNAseq']
            source = colidx['SmProtID']
            start = colidx['startCodonSite']
            stop = colidx['stopCodonSite']
            tags = [colidx['dataSource'], colidx['IsHighConfidence']]
            attributes = [colidx['transcriptID'], colidx['corresponding_EnsemblTranscriptIDs'],
                          colidx['corresponding_transcriptBiotypes'],
                          colidx['transcriptDNAseq'], colidx['corresponding_EnsemblGeneIDs'],
                          colidx['corresponding_geneBiotypes']]

        elif self.label == 'negative':
            seq = colidx['DNAseq']
            source = colidx['EnsemblTranscriptID']
            start = colidx['startCodonSite']
            stop = colidx['stopCodonSite']
            tags = []
            attributes = [colidx['transcriptBiotype'], colidx['transcriptDNAseq'],
                          colidx['EnsemblGeneID'], colidx['geneBiotype']]

        for line in infile:
            field = line.rstrip().split(',')
            entry = DataEntry(orf_id=field[seq],
                              orf_seq=field[source],
                              start_at=int(field[start]),
                              end_at=int(field[stop]),
                              class_label=self.label,
                              )
            self.append(entry)
            for t in tags:
                entry.tags.append(f'{col[t]}={field[t]}')
            for attr in attributes:
                entry.data.append(f'{col[attr]}={field[attr]}')

            self.append(entry)

        return len(self)

    def read_fasta(self, filename):
        """-----------------------------------------------------------------------------------------
        Read in sequences for prediction. Input is standard Fasta format, e.g.
        >SPROHSA001781
        ATGTATACGCTGCCTCGCCAGGCCACACCAGGTGTTCCTGCACAGCAGTCCCCAAGCATGTGA
        >SPROHSA001792
        ATGTGTGGTAACACCATGTCTGTGCCCCTGCTCACCGATGCTGCCACCGTGTCTGGAGCTGAGC

        All orfs with a start codon in DataSet.allowed_start_codons and a stop codon in
        DataSet.allowed_stop_codns are extracted. Not that RFs starting at a start codon but with
        no stop codon before the end of the sequence are NOT extracted (see ORF.py). The stop codon
        is not included in the extracted sequence.

        :param filename     string, sequence file in FastA format
        :return: int        number of orfs read
        -----------------------------------------------------------------------------------------"""
        infile = open(filename, 'r')
        if not infile:
            sys.stderr.write(f'DataSet:read_fasta() cannot open input file ({filename}')
            exit(2)

        sequences = []
        for line in infile:
            # read all sequences and store in dict indexed by sequence ID
            if line.startswith('>'):
                seqentry = {}
                sequences.append(seqentry)
                try:
                    id, doc = line.split(' ')
                except ValueError:
                    id = line
                    doc = ''
                seqentry['id'] = id.rstrip().replace('>', '')
                seqentry['doc'] = doc
                seqentry['seq'] = ''
            else:
                seqentry['seq'] += line.rstrip()

        # Convert to ORFS
        for seq in sequences:
            splitorf = ORF(seq['seq'], seq['id'], split=True)
            print(f"{seq['id']}\t{splitorf.pos}")
            for i, orf in enumerate(splitorf.pos):
                self.append(DataEntry(orf_id=f"{seq['id']}_{orf[0]}_{orf[1]}",
                                      orf_seq=seq['seq'][orf[0]:orf[1]],
                                      source=seq['id'],
                                      start_at=orf[0],
                                      end_at=orf[1],
                                      class_label=self.label
                                      ))

        return len(self)

    def orf_filter(self, label, constraints):
        """-----------------------------------------------------------------------------------------
        return a list of ORFs object with the ORFs of interest. Constraints can be used to filter
        the list. Constraints is a dictionary. The keys correspond to the column labels and the
        values are arrays that specify allowed values. For "DNAlength" the value indicates the
        minimum acceptable length (including stop codon).

        positive fields: SmProtID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,
        stopCodonSite,transcriptID,corresponding_EnsemblTranscriptIDs,
        corresponding_transcriptBiotypes,transcriptDNAseq,corresponding_EnsemblGeneIDs,
        corresponding_geneBiotypes,dataSource,IsHighConfidence

        negative fields: orfID,DNAseq,DNAlength,startCodon,stopCodon,startCodonSite,stopCodonSite,
        EnsemblTranscriptID,transcriptBiotype,transcriptDNAseq,EnsemblGeneID,geneBiotype

        :param self: DataSet
        :param label: string        should be "positive" or "negative"
        :param constraints: dict    keys=column labels, values = [] of allowed values
        :return: list of ORF        selected ORFs (see ORF.py)
        -----------------------------------------------------------------------------------------"""
        selected_orfs = []
        if label == "positive":
            for entry in self:
                orf = ORF(entry['DNAseq'], seqid=entry['SmProtID'])
                # DNAlength does not include stop codon
                orf.pos = [[0, int(entry['stopCodonSite']) - int(entry['startCodonSite']) + 3]]
                orf.offset = entry['startCodonSite']
                orf.label = 'positive'
                orf.tag = DataSet.to_unique_tag(entry['corresponding_transcriptBiotypes'])
                if entry['IsHighConfidence'] == 'Yes':
                    orf.tag.append('HighConfidence')

                if self.apply_constraints(entry, constraints):
                    selected_orfs.append(orf)

        elif label == "negative":
            negative_id_list = {}
            for entry in self:
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
                orf.tag = DataSet.to_unique_tag(entry['transcriptBiotype'])
                orf.split_orfs()

                if self.apply_constraints(entry, constraints):
                    selected_orfs.append(orf)

        else:
            sys.stderr.write(f'mipepid:DataSet:orf_filter() - unkown dataset label ({label})')

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


# ==================================================================================================
# End of class DataSet
# ==================================================================================================


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

    # fasta = DataSet('../datasets/negative_original.fasta_test.fa', 'unknown')

    # pos = DataSet(sys.argv[1], 'positive')
    # sys.stderr.write(f'Positive sequences read: {len(pos)}\n')
    # filtered_orfs = pos.orf_filter('positive', {'DNAlength': filter_len})
    # n_pos_filt = len(filtered_orfs)
    # sys.stderr.write(f'Positive sequences after filtering ({filter_len}): {n_pos_filt}\n')
    #
    neg = DataSet(sys.argv[2], 'negative')
    sys.stderr.write(f'Negative sequences read: {len(neg)}\n')
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

    # open output file and write column header
    columns = ['sORF_ID', 'sORF_seq', 'start_at', 'end_at', 'true_label', 'tags', 'classification', 'score',
               'probability']
    outfile = open(output_fname, 'w')
    outfile.write(','.join(columns))
    outfile.write('\n')
    outfile.close()

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
