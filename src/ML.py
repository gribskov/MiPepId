import sys
import numpy as np
import pickle
import pandas as pd
import warnings


# from sklearn.exceptions import InconsistentVersionWarning
# warnings.simplefilter("error", InconsistentVersionWarning)


class KmerFeaturization:
    """=============================================================================================
    Manage kmer features and conversion from sequence
    ============================================================================================="""

    def __init__(self, k):
        """-----------------------------------------------------------------------------------------
        k: length of kmer, typically 4
        -----------------------------------------------------------------------------------------"""
        self.k = k
        self.letters = ['A', 'T', 'C', 'G']
        # to convert to base 4 number, powers for the digits
        self.multiplyBy = 4 ** np.arange(k - 1, -1, -1)
        self.n = 4 ** k  # number of possible k-mers

    def seqlist_to_kmer(self, seqs, n_occur=False):
        """-----------------------------------------------------------------------------------------
        Convert a list of m DNA sequences to a 2-d array with shape (m, 4**k) for the 1-hot
        representation of the kmer features.
        
        :param seqs: list of string     DNA sequences to featurize
        :param n_occur: bool            if True the number of occurrences is saved
        :return: np.array               array of kmer features, shape=(m, 4**k)
        -----------------------------------------------------------------------------------------"""
        kmer_features = []
        for seq in seqs:
            this_kmer_feature = self.seq_to_kmer(seq.upper(), n_occur=n_occur)
            kmer_features.append(this_kmer_feature)

        kmer_features = np.array(kmer_features)

        return kmer_features

    def seq_to_kmer(self, seq, n_occur=False):
        """-----------------------------------------------------------------------------------------
        Given a DNA sequence, return the 1-hot representation of its kmer feature.
        
        :param seq: string      a single DNA sequence, upper case
        :param n_occur: bool    if True the number of occurrences is saved
        :return: 
        -----------------------------------------------------------------------------------------"""
        number_of_kmers = len(seq) - self.k + 1

        kmer_feature = np.zeros(self.n)

        for i in range(number_of_kmers):
            this_kmer = seq[i:(i + self.k)]
            this_numbering = self.kmer_to_onehot(this_kmer)
            kmer_feature[this_numbering] += 1

        if not n_occur:
            kmer_feature = kmer_feature / number_of_kmers

        return kmer_feature

    def kmer_to_onehot(self, kmer):
        """-----------------------------------------------------------------------------------------
        Given a k-mer, return its numbering (the 0-based position in 1-hot representation)
        
        :param kmer: 
        :return: 
       ----------------------------------------------------------------------------------------- """
        digits = []
        for letter in kmer:
            digits.append(self.letters.index(letter))

        digits = np.array(digits)
        numbering = (digits * self.multiplyBy).sum()

        return numbering


# ==================================================================================================
# End of class kmer_featurization
# ==================================================================================================

def predict(logr, x, threshold):
    """---------------------------------------------------------------------------------------------
    Predict whether a sequence is coding or non-coding
    
    :param logr: object         LogisticRegression 
    :param x:                   data for prediction (kmer features)
    :param threshold: float     probability for prediction, p>threshold is coding
    :return: int                y_pred, 0/1 indicating non-coding/coding
             float              y_pred_score, the prediction score value
             float              y_prob, the probability that an instance is in the assigned category
    ---------------------------------------------------------------------------------------------"""
    y_pred_score = logr.predict_proba(x)[:, 1]
    y_pred = (y_pred_score > threshold) + 0
    # the probability that an instance is in the assigned category
    y_prob = abs(1 - y_pred - y_pred_score)

    return y_pred, y_pred_score, y_prob


def load_model(model_fname='./src/model/model.pkl', no_model_warning=True):
    """---------------------------------------------------------------------------------------------
    Load the pickled LogisticRegression model. This is unstable as future version os sklearn may
    change the class making the pickled file unreadable. This what happened with the original 
    model.
    
    pickle.DEFAULT_PROTOCOL = 4 seems to work
    
    :param model_fname: string      path to a file with the model
    :param no_model_warning: bool   if True, warning about model version mismatch is surpressed
    :return: object                 logr
             float                  threshold for prediction
    ---------------------------------------------------------------------------------------------"""
    if no_model_warning:
        # disable sklearn warnings about model mismatch
        warnings.filterwarnings("ignore")

    f = open(model_fname, 'rb')
    logr = pickle.load(f)
    threshold = pickle.load(f)
    # model_v = logr.__getstate__()['_sklearn_version']
    # k = logr.__dict__.keys()
    f.close()

    return logr, threshold


def batch_predict(orfs, logr, threshold, output_fname, k=4):
    """---------------------------------------------------------------------------------------------
    Predict non-coding/coding on one batch of sequences and save to a dataframe. Prediction details
    sORF_ID                     ID of the ORF, if there are multiple ORFs in the original sequence
                                it will have suffices such as _ORF01, _ORF02 ...
    sORF_seq                    sequence of the extracted ORF
    transcript_DNA_sequence_ID  ID of the original sequence before extracting ORFs
    start_at                    beginning of the orf in the original sequence
    end_at                      end of the orf in the original sequence

    starting and ending positions are converted to 1-based for molecular biologists
    
    :param orfs:
    :param logr: object             LogisticRegression (from load_model)
    :param threshold: float         probability for prediction, p>threshold is coding
    :param output_fname: string     path to output file
    :param k: int                   kmer size
    :return: DataFrame              prediction details
    ---------------------------------------------------------------------------------------------"""
    class_dic = {
        1: 'coding',
        0: 'noncoding'
    }

    # start_at and end_at are both 1-based.
    columns = ['sORF_ID', 'sORF_seq', 'start_at', 'end_at', 'true_label', 'tags']
    df = pd.DataFrame(orfs, columns=columns)

    seqs = df.sORF_seq.tolist()
    obj = KmerFeaturization(k)
    kmer_features = obj.seqlist_to_kmer(seqs, n_occur=False)
    y_pred, _, y_prob = predict(logr, kmer_features, threshold)
    df['classification'] = [class_dic[x] for x in list(y_pred)]
    df['probability'] = y_prob
    df.to_csv(output_fname, header=True, index=False, mode='a')

    return df
