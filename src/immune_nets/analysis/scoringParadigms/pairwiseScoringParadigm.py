from itertools import combinations
import numpy as np

from immune_nets.analysis.statDistanceTypes import PairwiseStatDistance
from .scoringParadigm import ScoringParadigm

class PairwiseScoringParadigm(ScoringParadigm):
    def __init__(self, statDistance: PairwiseStatDistance):
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


