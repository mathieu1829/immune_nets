from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from src.creation.algorithms.common_methods import numerizeTCRSeq
import numpy as np
import math

class groupNegativeLevenshteinDistance:
    def tcr_dist(self,tcr):
        other = numerizeTCRSeq(tcr)
        ceil_all = np.vectorize(math.ceil) 
        return -ceil_all(squareform(pdist(other,metric="hamming"))*other.shape[1])

    def __str__(self):
        return "group_negative_levenshtein"

