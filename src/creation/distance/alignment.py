from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

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
        if ref != None:
            self_score = pairwise2.align.localds(ref, ref, self.__aligner.substitution_matrix, -10, -1)
        if x != None:
            alignments = pairwise2.align.localds(x, ref, self.__aligner.substitution_matrix, -10, -1)
            return 1 - (alignments[0].score / self_score) 
    def __str__(self):
        return self.__matrix

