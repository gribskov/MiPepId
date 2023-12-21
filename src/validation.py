"""=================================================================================================
MiPepid:validation.py
Original code from from QIHUA LIANG <qili00002@stud.uni-saarland.de>

mipepid.py output is csv with columns
sORF_ID,sORF_seq,start_at,end_at,true_label,tags,classification,score,probability


21 December 2023     gribskov
================================================================================================="""
import sys
import pandas as pd


def get_seq(i):
    """
    Does not seem needed
    :param i:
    :return:
    """
    a = input_file(i)
    file_out = './data/MiPepid/MiPepid_' + i + '_data.fa'
    print('output={}'.format(file_out))
    with open(file_out, 'w') as f_out:
        for id, seq in enumerate(a):
            f_out.write('>' + str(id) + '\n')
            f_out.write(seq + '\n')


def print_result():
    db = {'positive': pos, 'negative': neg}

    # print('\nLength >= {}'.format(length_cutoff))
    tf = {}

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


def stats(df):
    """---------------------------------------------------------------------------------------------

    :param df:
    :return: float, float     precision, recall
    ---------------------------------------------------------------------------------------------"""
    p = (df['true_label'] == 'positive').sum()
    n = (df['true_label'] == 'negative').sum()
    tp = (df['true_label'] == 'positive') & (df['classification'] == 'coding')
    fp = (df['true_label'] == 'negative') & (df['classification'] == 'coding')
    tp = tp.sum()
    fp = fp.sum()

    return p, n, tp / (tp + fp), tp / p


# ===================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    print('mipepid validation')
    try:
        threshold = float(sys.argv[2])
    except IndexError:
        threshold = 0.6
    lengths = [0, 110, 10]
    print(f'Probability threshold: {threshold}')
    print(f'Lengths: {lengths[0]} to {lengths[1]} by {lengths[2]}')

    # read the mipepid output
    infile_name = sys.argv[1]
    print(f'\nReading mipepid output from {infile_name}')
    df = pd.read_csv(infile_name, header=0, sep=',')
    df = df.astype({"start_at": "int", "end_at": "int"})
    print(f'{len(df)} read')

    # recall classification using desire P threshold
    pclass = ['noncoding', 'coding']
    df['classification'] = df['probability'].apply(lambda x: 'coding' if x > t else 'noncoding')

    print(f'\nP   len npos    nneg      prec  recall')
    for l in range(lengths[0], lengths[1], lengths[2]):
        # select by length
        len_ok = df['end_at'] - df['start_at'] + 1 > l
        df1 = df[len_ok]

        npos, nneg, precision, recall = stats(df1)
        print(f'{threshold}\t{l:3d}\t{npos:4d}\t{nneg:4d}\t{precision:.4f}\t{recall:.4f}')

exit(0)
