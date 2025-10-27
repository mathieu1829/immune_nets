from src.analysis.methods.graphStats import GraphStats
from scipy.spatial.distance import euclidean

class EuclideanStatDistance:
    def __init__(self):
        pass

    def stat_dist(self, a:GraphStats, b:GraphStats):
        return euclidean(a.toStatVector(),
                         b.toStatVector())
