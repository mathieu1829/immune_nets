from src.entities import GraphStats
from src.analysis.statDistanceTypes import PairwiseStatDistance
from scipy.spatial.distance import euclidean

class EuclideanStatDistance(PairwiseStatDistance):
    def __init__(self):
        pass

    def stat_dist(self, a:GraphStats, b:GraphStats) -> float:
        return euclidean(a.toStatVector(),
                         b.toStatVector())
