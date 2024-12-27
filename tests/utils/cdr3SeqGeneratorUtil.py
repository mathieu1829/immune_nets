import pandas as pd
import numpy as np
import random
from src.creation.distance.alignment import sequenceAligner

class cdr3SeqGeneratorUtil:
    def generate_random_cdr3(self,min_length = 11,max_length=16):
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        random_length = random.randint(min_length, max_length)
        return ''.join(random.choices(amino_acids, k=random_length))

    def mutateCDR3(self, cdr3seq, numMutations=1):
        seq = list(cdr3seq)
        for _ in range(numMutations):
            mutationIdx = random.randint(0,len(seq)-1)
            seq[mutationIdx] = random.choice("ACDEFGHIKLMNPQRSTVWY".replace(seq[mutationIdx],""))
        return "".join(seq)





def createRepertoireDataFrame(self, tcrA=None, tcrB=None, frequencies=None, n=None):
    if tcrA is None:
        if tcrB is None and n is None :
            raise ValueError("Length must be specified to generate repertoire")
        if tcrB is not None and n is not None :
            raise ValueError("n must not be set if tcrB repertoire set")

        if tcrB is not None :
            tcra = self.generateRelatedSubset







if __name__ == "__main__":
    testGen = cdr3SeqGeneratorUtil()
    randomCDR3 = testGen.generate_random_cdr3()
    print(f"random cdr3: {randomCDR3}")
    print("Mutated random cdr3:")
    for _ in range(10):
        print(testGen.mutateCDR3(randomCDR3))
    
