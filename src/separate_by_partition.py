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
    threshold = 0.6
    data = 'datasets/'
    # partition[test/train][label] = []
    partition = {'train': [[], []],
                 'test':  [[], []]}
    test = 'test'
    train = 'train'

    # read training sequence IDs from training_feature.csv
    # column 256 and 257 are ORF_ID and label (1/0) respectively
    # training = []
    train_in = open(data + 'train_feature.csv', 'r')
    header = train_in.readline()
    for line in train_in:
        field = line.rstrip().split(',')
        orf = field[256]
        label = field[257]
        print(f'{label}\t{orf}')
        partition[train][int(label)].append(orf)

    train_in.close()
    
    # read test sequence IDs from test_feature.csv
    # testing = []
    test_in = open(data + 'test_feature.csv', 'r')
    header = test_in.readline()
    for line in test_in:
        field = line.rstrip().split(',')
        orf = field[256]
        label = field[257]
        print(f'{label}\t{orf}')
        partition[test][int(label)].append(orf)

    test_in.close()

    print(f'\npartition\tpositive\tnegative')
    print(f'training\t{len(partition[train][1]):8d}\t{len(partition[train][0]):8d}')
    print(f'test    \t{len(partition[test][1]):8d}\t{len(partition[test][0]):8d}')
    print(f'total   \t{len(partition[test][1])+len(partition[train][1]):8d}\t{len(partition[test][0])+len(partition[train][0]):8d}')

    # make an index of all sequences
    idx = {}
    for part in partition:
        for posneg in (0,1):
            for seq in partition[part][posneg]:
                idx[seq] = [part, posneg]


    # partition positive results into training/test partition
    # sORF_ID, sORF_seq, transcript_DNA_sequence_ID, start_at, end_at, classification, probability
    pos = open('positive_results.csv', 'r')
    header = pos.readline()
    missing = 0
    for line in pos:
        field = line.rstrip().split(',')
        orig_orf = orf
        orf = field[0].replace('ORF', 'ORF0')
        classid = field[5]
        p = float(field[6])
        call = 'noncoding'
        if p > threshold:
            call = 'coding'
        if orf in idx:
            idx[orf] += [orig_orf, classid, call, p]
        else:
            cut = orf.index('_')
            orig_orf = orf
            orf = orf[:cut]
            if orf in idx:
                idx[orf] += [orig_orf, classid, call, p]
            else:
                missing += 1
                print(f'{missing}\t{orf} is missingclass={classid}, p={p}')

    # partition negative results into training/test and by biotype

    n = 0
    for id in idx:
        n += 1
        print(f'{n}\t{id}\t{idx[id]}')
    exit(0)
