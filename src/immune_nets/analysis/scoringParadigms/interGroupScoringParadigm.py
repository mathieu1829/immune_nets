from itertools import combinations
import numpy as np

from immune_nets.analysis.statDistanceTypes import InterGroupStatDistance
from .scoringParadigm import ScoringParadigm

class InterGroupScoringParadigm(ScoringParadigm):
    def __init__(self, statDistance: InterGroupStatDistance):
        self.statDistance = statDistance
    
    def compute_score(self, groupedStats):
        groups = [group for group in groupedStats]
        inter_group_distances = []
        for combo in combinations(groups,2):
            group_a = groupedStats[combo[0]]
            group_b = groupedStats[combo[1]]
            score = self.statDistance.group_dist(group_a, group_b)
            inter_group_distances.append(score)
        inter_group_distances = np.array(inter_group_distances)
        return inter_group_distances.mean() - np.std(inter_group_distances) 



