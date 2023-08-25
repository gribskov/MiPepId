"""=================================================================================================
write out the fasta sequences for the test and training partitions based on the IDs in the \test_feature.csv and
training_feature.csv

Michael Gribskov     25 August 2023
================================================================================================="""
import sys

# --------------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    data = '../datasets/'

    # read training sequence IDs from training_feature.csv
    # column 257 and 258 are ORF_ID and label (1/0) respectively
    train = []
    train_in = open( data + 'train_feature.csv', 'r' )
    header = train_in.readline()
    for line in train_in:
        field = line.rstrip().split(',')
        orf = field[257]
        label = field[258]
        print(f'{label}\t{orf}')


    # read test sequence IDs from test_feature.csv
    test = []


    # partition positive results into training/test partition

    # partition negative results into training/text and by biotype

    exit(0)
