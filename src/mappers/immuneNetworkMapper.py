import pandas as pd
import pickle

from src.creation.immuneNetwork import ImmuneNetwork

class ImmuneNetworkMapper:
    @staticmethod
    def toPickle(network: ImmuneNetwork, path):
        with open(path,"wb") as f:
            pickle.dump(network,f)
    
