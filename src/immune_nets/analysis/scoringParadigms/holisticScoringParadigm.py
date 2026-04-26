from itertools import combinations
import numpy as np

from immune_nets.analysis.statDistanceTypes import HolisticStatDistance 
from .scoringParadigm import ScoringParadigm

class HolisticScoringParadigm(ScoringParadigm):
    """
    Computes a score over all groups in a single global computation.

    This paradigm is used when the distance is defined directly on the full
    grouped dataset, rather than being decomposed into pairwise or inter-group
    comparisons.
    """
    def __init__(self, statDistance: HolisticStatDistance):
        """
        :param statDistance: Computes a scoredirectly from the full groupedStats 
                             structure without decomposing into pairwise or 
                             inter-group comparisons.
        """

        self.statDistance = statDistance

    def compute_score(self, groupedStats):
        return self.statDistance.compute_holistic_score(groupedStats)
