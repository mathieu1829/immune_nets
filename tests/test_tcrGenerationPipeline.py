import unittest
from tests.utils.tcrGenerationPipeline import *
from tests.utils.frequencyGenerator import *
from unittest.mock import MagicMock
from src.creation.distance.alignment import sequenceAligner

class TestTcrGenerationPipeline(unittest.TestCase):

    def test_error_handling(self):
        tcrPipeline = tcrGenerationPipeline(1)
        freqGen = frequencyGenerator().generateLowFrequencies
        repRange = MagicMock()
        repRange._length = 1

        def returnMockLength():
            return repRange._length

        repRange.__len__.side_effect = returnMockLength

        # Testing if incorrect ranges will be accepted
        try:
            tcrPipeline._validateGeneratorRepertoireRange(
                    n=1,
                    idx = 1,
                    rangeStart=1
                    )
            self.fail("It should not be able to use both indexes and ranges in repertoire range")
        except:
            pass

        try:
            tcrPipeline._validateGeneratorRepertoireRange(
                    n=1,
                    rangeStart=1,
                    )
            self.fail("Both threshold for the range must be provided if thresholds are used")
        except:
            pass

        # Testing alternate tcr parameters

        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=freqGen,
                    tcrA = "aaaa",
                    repertoireRange=repRange
                    )
            self.fail("Alternative tcrA sequence should not be able to be passed as other types than method or array")
        except:
            pass

        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=freqGen,
                    tcrB = "aaaa",
                    repertoireRange=repRange
                    )
            self.fail("Alternative tcrB sequence should not be able to be passed as other types than method or array")
        except:
            pass

        # Testing frequencies 
        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=freqGen
                    )
            self.fail("Repertoire Range parameter should be present but it is not")
        except:
            pass

        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=" ",
                    repertoireRange=repRange
                    )
            self.fail("""Frequencies should be one of:
                                * method
                                * array of size n
                                * list of size of repertoire range of numpy arrays of size n, 
                                * list of size of repertoire range of methods
                            """)
        except:
            pass

        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=np.ones((12,2)),
                    repertoireRange=repRange
                    )
            self.fail("If frequencies are passed as an array, array should be one-dimensional")
        except:
            pass

        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=np.ones(12),
                    repertoireRange=repRange
                    )
            self.fail("If frequencies are passed as an array length of an array should be equal n")
        except:
            pass

        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=[np.ones(1), np.ones([1])],
                    repertoireRange=repRange
                    )
            self.fail("If frequencies are passed as an list, the size of the list should be equal length of repertoire range")
        except:
            pass
        repRange._length = 2

        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=[np.array([1]),freqGen],
                    repertoireRange=repRange
                    )
            self.fail("If frequencies are passed as a list, all elements inside must be of the same type")
        except:
            pass

        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=[np.ones(1),np.ones((2,2))],
                    repertoireRange=repRange
                    )
            self.fail("If frequencies are passed as a list or arrays, all arrays should be one dimensional")
        except:
            pass

        try:
            tcrPipeline._validateGeneratorParameters(
                    n=1,
                    frequencies=[np.ones(1),np.ones(2)],
                    repertoireRange=repRange
                    )
            self.fail("If frequencies are passed as a list or arrays, all arrays should have length equal to n")
        except:
            pass
    
    def test_randomRepertoireGeneration(self):
        tcrPipeline = tcrGenerationPipeline(1)
        freqGen = frequencyGenerator()
        rogueTcr = np.array(["AAAAAAAAAAA"])
        # print(tcrPipeline.generateRandomRepertoire(n=5,frequencies=freqGen.generateLowFrequencies)
        #       .generateRandomRepertoire(n=3,frequencies=freqGen.generateVeryHighFrequencies)
        #       .generateRandomRepertoire(n=1,frequencies=freqGen.generateSimulatedNaturalFrequencies,tcrA=rogueTcr)
        #       .generateRandomRepertoire(n=1,frequencies=freqGen.generateSimulatedNaturalFrequencies,tcrB=rogueTcr)
        #       .generateRandomRepertoire(n=1,frequencies=freqGen.generateSimulatedNaturalFrequencies,tcrA=rogueTcr,cloneTcrInclusive=True)
        #       .toDataFrame()[0])

        manyTcrPipeline = tcrGenerationPipeline(3)
        # dfs = manyTcrPipeline.generateRandomRepertoire(n=5,frequencies=freqGen.generateLowFrequencies) \
        #         .generateRandomRepertoire(n=2,frequencies=freqGen.generateLowFrequencies,idx=0) \
        #         .generateRandomRepertoire(n=2,frequencies=freqGen.generateVeryHighFrequencies,rangeStart=1,rangeEnd=2) \
        #         .generateRandomRepertoire(n=1,frequencies=freqGen.generateMediumFrequnencies, cloneTcrA=True) \
        #         .generateRandomRepertoire(n=1,frequencies=freqGen.generateMediumFrequnencies, cloneTcrB=True) \
        #         .generateRandomRepertoire(n=1,frequencies=freqGen.generateMediumFrequnencies, cloneTcrB=True, cloneTcrInclusive=False) \
        #         .generateRandomRepertoire(n=1,frequencies=freqGen.generateMediumFrequnencies, cloneTcrB=True, cloneTcrA=True) \
        #         .toDataFrame()

        manyTcrPipeline.generateRelatedRepertoire(n=20,frequencies=freqGen.generateHighFrequencies,distance=0.15, distanceFun=sequenceAligner("BLOSUM62").tcr_dist, split=True) \
                .generateRelatedRepertoire(n=5,frequencies=freqGen.generateSimulatedNaturalFrequencies,distance=0.15, distanceFun=sequenceAligner("BLOSUM62").tcr_dist) \
                .generateRandomRepertoire(n=50,frequencies=freqGen.generateSimulatedNaturalFrequencies) \
                .generateRandomRepertoire(n=2,frequencies=freqGen.generateHighFrequencies,idx=0) \
                .generateRandomRepertoire(n=2,frequencies=freqGen.generateHighFrequencies,idx=2) \
                .generateRandomRepertoire(n=2,frequencies=freqGen.generateVeryHighFrequencies,idx=1) 

        for _ in range(5):
            manyTcrPipeline.generateRelatedRepertoire(n=15,frequencies=freqGen.generateSimulatedNaturalFrequencies,distance=0.15, distanceFun=sequenceAligner("BLOSUM62").tcr_dist,split=True) 
        for _ in range(10):
            manyTcrPipeline.generateRelatedRepertoire(n=10,frequencies=freqGen.generateSimulatedNaturalFrequencies,distance=0.15, distanceFun=sequenceAligner("BLOSUM62").tcr_dist,split=True) 
        for _ in range(5):
            manyTcrPipeline.generateRelatedRepertoire(n=5,frequencies=freqGen.generateSimulatedNaturalFrequencies,distance=0.15, distanceFun=sequenceAligner("BLOSUM62").tcr_dist,split=True) 
        for _ in range(10):
            manyTcrPipeline.generateRelatedRepertoire(n=3,frequencies=freqGen.generateSimulatedNaturalFrequencies,distance=0.15, distanceFun=sequenceAligner("BLOSUM62").tcr_dist,split=True) 

        dfs = manyTcrPipeline.toDataFrame()


        # for df in dfs:
        #     print(df)
        # for i,df in enumerate(dfs):
        #     df.to_csv(f"publicTest{i}.csv")
        #
        

    

