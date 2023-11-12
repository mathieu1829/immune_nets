from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

class sequenceAligner:
    def __init__(self, matrix):
        self.aligner = PairwiseAligner()
        self.aligner.substitution_matrix = substitution_matrices.load(matrix)

    def tcr_alig(self,x, ref):
        if x != None:
            alignments = pairwise2.align.localds(x, ref, self.aligner.substitution_matrix, -10, -1)
            return alignments[0].score
