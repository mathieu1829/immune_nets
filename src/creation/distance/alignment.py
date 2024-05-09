from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import numpy as np

class sequenceAligner:
    def setMatrix(self,matrix):
        self.__matrix = matrix
        self.__aligner = PairwiseAligner()
        self.__aligner.substitution_matrix = substitution_matrices.load(matrix)
        self.__aligner.open_gap_score = -10;
        self.__aligner.extend_gap_score = -1;

    def getMatrix(self):
        return self.__matrix

    def __init__(self, matrix,switch_self_score=False):
        self.setMatrix(matrix)
        self.switch_self_score = switch_self_score

        
        

    def tcr_dist(self,x, ref):

        assert x != None, f"x is None, can't calculate distace"
        assert ref != None, f"ref is None, can't calculate distace"
        
        self_score_x = self.__aligner.align(x,x)
        self_score_ref = self.__aligner.align(ref,ref)

        self_score = max(self_score_x.score, self_score_ref.score)
        alignments = self.__aligner.align(x, ref)

        return 1 - (float(alignments.score) / float(self_score)) 
    def __str__(self):
        return self.__matrix

