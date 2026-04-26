from itertools import combinations
import numpy as np

from immune_nets.analysis.statDistanceTypes import PairwiseStatDistance
from .scoringParadigm import ScoringParadigm

class PairwiseScoringParadigm(ScoringParadigm):
    """
    Computes a score by performing pairwise comparisons between samples of groups.

    For each pair of groups, all pairwise distances between their samples are computed.
    These values are aggregated into a single group-to-group distance. The final score
    is the mean of all inter-group distances penalized by their standard deviation.
    """

    def __init__(self, statDistance: PairwiseStatDistance):
        """
        :param statDistance: Computes a distance between two individual GraphStats
                             instances (samples).
        """
        self.statDistance = statDistance

    def compute_score(self, groupedStats):
        groups = [group for group in groupedStats]
        inter_group_distances = []
        for combo in combinations(groups,2):
            group_a = groupedStats[combo[0]]
            group_b = groupedStats[combo[1]]
            inter_group_distance = np.array([self.statDistance.stat_dist(a,b) for a in group_a for b in group_b ])
            inter_group_distances.append(inter_group_distance.mean())
        inter_group_distances = np.array(inter_group_distances)
        return inter_group_distances.mean() - np.std(inter_group_distances) 


