"""=================================================================================================
write out the fasta sequences for the test and training partitions based on the IDs in the \test_feature.csv and
training_feature.csv

Michael Gribskov     25 August 2023
================================================================================================="""
import sys

def read_features(fname, ptype, partition):
    '''---------------------------------------------------------------------------------------------
    read in entire feature list and store in partition
     partition = {'train': [[], []],
                 'test':  [[], []]}
    :param fname:
    :param
    :param partition:
    :return:
    ---------------------------------------------------------------------------------------------'''

# --------------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    data = 'datasets/'
    # partition[test/train][label] = []
    partition = {'train': [[], []],
                 'test':  [[], []]}
    test = 'test'
    train = 'train'

    # read training sequence IDs from training_feature.csv
    # column 256 and 257 are ORF_ID and label (1/0) respectively
    train = []
    train_in = open(data + 'train_feature.csv', 'r')
    header = train_in.readline()
    for line in train_in:
        field = line.rstrip().split(',')
        orf = field[256]
        label = field[257]
        print(f'{label}\t{orf}')
        partition[test][int(label)].append(orf)

    # read test sequence IDs from test_feature.csv
    test = []

    # partition positive results into training/test partition

    # partition negative results into training/text and by biotype

    exit(0)
