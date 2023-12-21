"""=================================================================================================
MiPepid:validation.py
Original code from from QIHUA LIANG <qili00002@stud.uni-saarland.de>


21 December 2023     gribskov
================================================================================================="""

import sys
import pandas as pd


def input_file(i):
    df = pd.read_csv('./data/MiPepid/datasets/' + i + '_original_data.csv', header=0, sep=',')
    seq = df['DNAseq']
    seq.rename('RNAseq')
    return seq


def get_seq(i):
    a = input_file(i)
    file_out = './data/MiPepid/MiPepid_' + i + '_data.fa'
    print('output={}'.format(file_out))
    with open(file_out, 'w') as f_out:
        for id, seq in enumerate(a):
            f_out.write('>' + str(id) + '\n')
            f_out.write(seq + '\n')


def print_result():
    threshold = 0.75
    pos = sys.argv[1]
    neg = sys.argv[2]
    print('positive data:{}'.format(pos));
    print('negative data:{}'.format(neg));
    db = {'positive': pos, 'negative': neg}
    for length_cutoff in range(0, 110, 10):
        # print('\nLength >= {}'.format(length_cutoff))
        tf = {}
        print('\nlength:{}\tthreshold:{}'.format(length_cutoff, threshold))
        for i in db:
            df = pd.read_csv(db[i])
            select = (df['end_at'] - df['start_at'] + 1) >= length_cutoff
            ii = i[:3]
            # print('{} -> {}'.format(i, ii))
            # print(select)
            ##print (df['sORF_ID'])
            # print (df['sORF_ID'].str.endswith('_ORF1'))
            # df1 = df.loc[df['sORF_ID'].str.endswith('_ORF1'), ['transcript_DNA_sequence_ID', 'classification']]
            df1 = df.loc[
                select, ['transcript_DNA_sequence_ID', 'classification', 'probability']]
            # df1['alt'] = df1['probability'] >= threshold
            df1['alt'] = df1['probability'].apply(
                lambda x: 'altcoding' if x > threshold else 'altnoncoding')
            # f1['A'] = df1['A'].apply(lambda x: [y if y <= 9 else 11 for y in x])
            pts = df1.alt.value_counts()
            print('    {:8s}{:8d}{:8d}'.format(ii, pts[0], pts[1]))
            # print(df1.head())
            cts = df1.classification.value_counts()
            # tf[ii] = [cts[0], cts[1]]
            tf[ii] = [pts[0], pts[1]]
            # print(cts)

        n = 'neg'
        p = 'pos'
        print('\t\tPos\tNeg\tAll')
        rp = tf[p][0] / (tf[p][0] + tf[p][1])
        rn = tf[n][1] / (tf[n][0] + tf[n][1])
        ra = (tf[p][0] + tf[n][1]) / (tf[p][0] + tf[p][1] + tf[n][0] + tf[n][1])
        pp = tf[p][0] / (tf[p][0] + tf[n][0])
        pn = tf[n][1] / (tf[p][1] + tf[n][1])
        pa = (tf[p][0] + tf[n][1]) / (tf[p][0] + tf[p][1] + tf[n][0] + tf[n][1])
        print('recall\t\t{:5.3f}\t{:5.3f}\t{:5.3f}'.format(rp, rn, ra))
        print('precision\t{:5.3f}\t{:5.3f}\t{:5.3f}'.format(pp, pn, pa))


# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    # 1. use get_seq() to extract DNA sequences
    # get_seq('positive')
    # get_seq('negative')
    # 2. run mipepid program in command line
    # python3 ./src/mipepid.py  ../data/MiPepid/MiPepid_positive_data.fa ../data/MiPepidOutput/MiPepid_pos_data1.csv
    # python3 ./src/mipepid.py  ../data/MiPepid/MiPepid_negative_data.fa ../data/MiPepidOutput/MiPepid_neg_data1.csv
    # 3. summarize the classification results
    print_result()
    ############# RESULT ###############
    # MiPepid_pos
    # coding       3907
    # noncoding      80
    # Name: classification, dtype: int64
    # MiPepid_neg
    # coding       2135
    # noncoding     801
    # Name: classification, dtype: int64
    ####################################

exit(0)
