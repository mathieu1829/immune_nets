import pickle
import pandas as pd

from src.entities import ImmuneNetwork

class ImmuneNetworkFactory:
    @staticmethod
    def fromPickle(path) -> ImmuneNetwork:
        with open(path, "rb") as f:
            immuneNet = pickle.load(f)
        return immuneNet

    @staticmethod
    def createEmpty():
        network = ImmuneNetwork(graph=pd.DataFrame({"r1":[],"r2":[]}),
                                proportions=[],
                                method="dummyMethod",
                                distanceFun="dummyDistanceFun",
                                sampleSize=0,
                                sampleId=None
                                )
        return network
