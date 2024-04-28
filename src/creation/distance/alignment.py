from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import numpy as np

class sequenceAligner:
    def setMatrix(self,matrix):
        self.__matrix = matrix
        self.__aligner = PairwiseAligner()
        self.__aligner.substitution_matrix = substitution_matrices.load(matrix)

    def getMatrix(self):
        return self.__matrix

    def __init__(self, matrix,switch_self_score=False):
        self.setMatrix(matrix)
        self.switch_self_score = switch_self_score

        
        

    def tcr_dist(self,x, ref):

        assert x != None, f"x is None, can't calculate distace"
        assert ref != None, f"ref is None, can't calculate distace"
        
        self_score_x = pairwise2.align.localds(x, x, self.__aligner.substitution_matrix, -10, -1)
        self_score_ref = pairwise2.align.localds(ref, ref, self.__aligner.substitution_matrix, -10, -1)
        self_score = max(self_score_x[0].score, self_score_ref[0].score)
        alignments = pairwise2.align.localds(x, ref, self.__aligner.substitution_matrix, -10, -1)

        return 1 - (float(alignments[0].score) / float(self_score)) 
    def __str__(self):
        return self.__matrix

