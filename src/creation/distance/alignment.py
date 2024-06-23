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

    def __init__(self, matrix,switch_self_score=False, group = False, negative = False):
        self.setMatrix(matrix)
        self.switch_self_score = switch_self_score
        self.group = group
        self.negative = negative

        
        

    def tcr_dist(self,x, ref = None):
        if x is None :
            raise ValueError('\"x\" parameter must not be none')

        if ref is not None:
            if self.group == True:
                raise ValueError('\"group\" must not be True, when using method in single comparison mode')
            self_score_x = self.__aligner.align(x,x)
            self_score_ref = self.__aligner.align(ref,ref)

            self_score = max(self_score_x.score, self_score_ref.score)
            alignments = self.__aligner.align(x, ref)

            res = 1 - (float(alignments.score) / float(self_score)) 
        else:
            if self.group == False:
                raise ValueError('\"group\" must be set to true when using method in group mode')
            tcr = x
            res = pdist(tcr,metrix = lambda u,v : "insert something here" )

        return res

    def __str__(self):
        return self.__matrix

