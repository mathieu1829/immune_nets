from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from scipy.spatial.distance import hamming
from scipy.spatial.distance import pdist
import numpy as np

def split_tcr_column(x, subunit):
    s1 = x.split(';')
    for sx in s1:
        s = sx.split(':')
        if s[0] == subunit:
            return s[1]

def equalizeTCRSeq(aa_tcr):
    longest = len(max(aa_tcr.tolist(),key = len))
    ljustToLongest = lambda x : x.ljust( longest)
    ljustToLongestVector = np.vectorize(ljustToLongest)
    return np.array([ list(i) for i in ljustToLongestVector(aa_tcr)])

def makeNumeric(aa_tcr):
    toInt = np.vectorize(ord)
    return toInt(aa_tcr)

def numerizeTCRSeq(aa_tcr):
    return makeNumeric(equalizeTCRSeq(aa_tcr))

if __name__ == "__main__":
    test = np.array(['aa','a','aaa'])

    print(equalizeTCRSeq(test))
    print(numerizeTCRSeq(test))

    # test = equalizeTCRSeq(test)

    # print(test)
    # print(hamming(test[0,],test[1,]))
    # print(pdist(test,metric='hamming'))

