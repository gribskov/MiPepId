"""=================================================================================================
Read the results file and calculate the number of positive orfs
result file format is
    sORF_ID,sORF_seq,transcript_DNA_sequence_ID,start_at,end_at,classification,probability
one row for each ORF

Michael Gribskov     23 August 2023
================================================================================================="""
import sys

def maxp_result(result):
    '''---------------------------------------------------------------------------------------------
    classify per sequence using the highest probability of any orf being coding

    :param result:
    :return:
    ---------------------------------------------------------------------------------------------'''
    header = result.readline()

    name_old = ''
    orf_n = 0
    p_max = 0
    seq_n = 0
    threshold = 0.6
    pos_n = 0
    neg_n = 0
    for line in result:
        (sorf_id, seq, t_id, start, end, call, p) = line.rstrip().split(',')
        p = float(p)

        if not name_old:
            # first sequence
            name_old = t_id

        if t_id == name_old:
            p_max = max( p_max, p)
        else:
            seq_n += 1
            print(f'{seq_n}\t{name_old}\t{p_max}')
            if p_max > threshold:
                pos_n += 1
            else:
                neg_n += 1
            name_old = t_id
            p_max = p

    print(f'{seq_n}\t{name_old}\t{p_max}')
    print(f'pos:{pos_n}\nneg:{neg_n}')

    return
# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    result = open(sys.argv[1], 'r')
    header = result.readline()

    orf_n = 0

    pos = 'pos'
    neg = 'neg'
    confusion = {'coding':{'n':0, pos:0, neg:0},
                 'noncoding':{'n':0, pos:0, neg:0}}


    orf_n = 0
    threshold = 0.6

    for line in result:
        (sorf_id, seq, t_id, start, end, call, p) = line.rstrip().split(',')
        p = float(p)

        confusion[call]['n'] += 1
        if p > threshold:
            confusion[call][pos] += 1
        else:
           confusion[call][neg] += 1

    for call in ('coding', 'noncoding'):
        print(f'call is {call}')
        print(f"{confusion[call]['n']}\t{confusion[call][pos]}\t{confusion[call][neg]}")



    exit(0)
