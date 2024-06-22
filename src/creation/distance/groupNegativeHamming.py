from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from src.creation.algorithms.common_methods import numerizeTCRSeq

class groupNegativeHammingDistance:
    def tcr_dist(self,tcr):
        tcr = numerizeTCRSeq(tcr)
        return -squareform(pdist(tcr,metric="hamming"))
        

    def __str__(self):
        return "group_hamming"


