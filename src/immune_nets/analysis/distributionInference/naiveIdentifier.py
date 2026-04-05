import numpy as np

from .distributionIdentifier import DistributionIdentifier

class NaiveIdentifier(DistributionIdentifier):
    def is_uniform_pmf(self, probs, tolerance=0.05):
        n = len(probs)
        target = np.ones(n) / n
        return np.all(np.abs(probs - target) < tolerance)

    def is_powerlaw_pmf(self, probs):
        if probs.argmax() == 0 and probs[0] >= 0.85:
            return True
        isDescending = np.all(np.diff(probs) <= 0)
        if probs.argmax() == 0 and isDescending and probs[0] > (1/len(probs)):
            return True 
        return False

    def is_normal_pmf(self, probs):
        maxIdx = probs.argmax()
        middle = []
        middleIdx = len(probs) // 2
        if len(probs) % 2:
            middle.append(middleIdx+1)
        else:
            middle.append(middleIdx)
            middle.append(middleIdx+1)

        leftSideAscending = np.all(np.diff(probs[:maxIdx]) >= 0)
        rightSideDescending = np.all(np.diff(probs[maxIdx:]) <= 0)
        if maxIdx in middle and leftSideAscending and rightSideDescending:
            return True   
        return False




    def identify_distribution(self, statDict):
        values = np.array(list(statDict.keys()))
        probs  = np.array(list(statDict.values()))
  
        if len(values) == 1:
            return "uniform"
        if self.is_normal_pmf(probs):
            return "normal"
        if self.is_powerlaw_pmf(probs):
            return "powerlaw"
        if self.is_uniform_pmf(probs):
            return "uniform"
        return "gamma"
