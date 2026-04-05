import pandas as pd
import pickle

from immune_nets.entities import ImmuneNetwork

class ImmuneNetworkMapper:
    @staticmethod
    def toPickle(network: ImmuneNetwork, path):
        with open(path,"wb") as f:
            pickle.dump(network,f)
    
