"""=================================================================================================
Using the assingment to test and training sets in datasets/test_feature.csv and train_feature.csv
produce positive and negative test and training fasta files. The sequence ID and the label (0/1)
are the last two columns in the csv file. Extract them with

cut -d ","  -f 257,258 datasets/train_feature.csv | grep ",0" |cut -d "," -f 1 > train.neg.id.txt

Michael Gribskov     11 December 2023
================================================================================================="""
import sys

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    datafiles = {'trainpos': '../datasets/train.pos.id.txt',
                 'trainneg': '../datasets/train.neg.id.txt',
                 'testpos':  '../datasets/test.pos.id.txt',
                 'testneg':  '../datasets/test.neg.id.txt',
                 }
    csvfiles = {'positive': '../datasets/positive_original.fasta',
                'negative': '../datasets/negative_original.fasta'}

    # open and read IDs
    data = {}
    label = 'pos'
    for f in datafiles:
        infile = open(datafiles[f], 'r')
        data[f] = []
        n = 0
        for line in infile:
            n += 1
            data[f].append(line.rstrip())

        print(f'{f}:{n}')
        infile.close()

    found = {d: 0 for d in data}

    label = 'pos'
    for c in csvfiles:
        infile = open(csvfiles[c], 'r')

        train_id = 'train' + label
        test_id = 'test' + label
        out = {train_id: open(csvfiles[c] + '_train.fa', 'w'),
               test_id:  open(csvfiles[c] + '_test.fa', 'w')}

        match_id = 'unknown'
        unmatched = 0
        matched = 0
        for line in infile:
            match_id = ''
            id = line.rstrip()[1:]
            for d in (train_id, test_id):
                if id in data[d]:
                    found[d] += 1
                    seq = infile.readline()
                    out[d].write(f'>{id}\n{seq}')
                    match_id = d
                    matched += 1
                    break

                if match_id == 'unknown':
                    unmatched += 1

            # print(f'{id} -> {match_id}')

        print(f'matched:{matched}')
        print(f'unmatched:{unmatched}')
        label = 'neg'

    exit(0)
