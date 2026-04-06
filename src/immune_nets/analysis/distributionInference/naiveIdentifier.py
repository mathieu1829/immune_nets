import numpy as np

from .distributionIdentifier import DistributionIdentifier

class NaiveIdentifier(DistributionIdentifier):
    """
    Naive distribution identifier based on simple heuristics.

    Classifies a probability mass function (PMF) by applying a sequence of
    rule-based checks:
    - If the PMF contains a single value, it is classified as uniform.
    - Otherwise, it checks for normality, then power-law behavior, then uniformity.
    - If none of the above match, the distribution is classified as gamma.
    """

    def is_uniform_pmf(self, probs, tolerance=0.05):
        """
        Check whether the PMF resembles a uniform distribution.

        Compares probabilities to a uniform distribution within a given tolerance.

        :param probs: list of probabilities in PMF
        :param tolerance: maximum tolerated divergence from uniform distribution
        """

        n = len(probs)
        target = np.ones(n) / n
        return np.all(np.abs(probs - target) < tolerance)

    def is_powerlaw_pmf(self, probs):
        """
        Check whether the PMF resembles a power-law distribution.

        Heuristic:
        - The first value has the highest probability, and
        - Probabilities show a generally decreasing trend.

        :param probs: list of probabilities in PMF
        """

        if probs.argmax() == 0 and probs[0] >= 0.85:
            return True
        isDescending = np.all(np.diff(probs) <= 0)
        if probs.argmax() == 0 and isDescending and probs[0] > (1/len(probs)):
            return True 
        return False

    def is_normal_pmf(self, probs):
        """
        Check whether the PMF resembles a normal distribution.

        Heuristic:
        - The highest probability occurs near the center, and
        - Probabilities increase up to the peak and decrease afterward.

        :param probs: list of probabilities in PMF
        """

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
