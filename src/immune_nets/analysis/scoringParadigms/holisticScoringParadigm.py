from itertools import combinations
import numpy as np

from immune_nets.analysis.statDistanceTypes import HolisticStatDistance 
from .scoringParadigm import ScoringParadigm

class HolisticScoringParadigm(ScoringParadigm):
    def __init__(self, statDistance: HolisticStatDistance):
        self.statDistance = statDistance

    def compute_score(self, groupedStats):
        return self.statDistance.compute_holistic_score(groupedStats)
