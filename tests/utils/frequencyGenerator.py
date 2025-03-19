import numpy as np

class frequencyGenerator:
    def generateLowFrequencies(self, n):
        options = np.arange(1,4)
        probabilities = [0.80,0.15, 0.05]
        return np.random.choice(a = options, p=probabilities, size=n)

    def generateMediumFrequnencies(self, n):
        options = np.arange(3,11)
        return np.random.choice(a = options, size=n)

    def generateHighFrequencies(self, n):
        options = np.arange(10,101)
        return np.random.choice(a = options, size=n)

    def generateVeryHighFrequencies(self, n):
        options = np.arange(101,1001)
        return np.random.choice(a = options, size=n)

    def generateSimulatedNaturalFrequencies(self, n):
        funList = np.array([self.generateLowFrequencies,self.generateMediumFrequnencies, self.generateHighFrequencies, self.generateVeryHighFrequencies])
        probabilities = [0.80, 0.12, 0.05, 0.03]
        return np.array([ np.random.choice(a=funList,p=probabilities)(1)[0] for _ in range(n) ])
