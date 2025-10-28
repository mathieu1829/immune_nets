import pickle

from src.creation.immuneNetwork import ImmuneNetwork

class ImmuneNetworkFactory:
    @staticmethod
    def fromCustomObject(obj) -> ImmuneNetwork:
        return ImmuneNetwork(graph=obj.graph,
                             method=obj.method,
                             distanceFun=obj.distanceFun,
                             sampleSize=obj.sampleSize,
                             threshold=obj.threshold,
                             threshold_alpha=obj.threshold_alpha,
                             threshold_beta=obj.threshold_beta,
                             name=obj.name,
                             sampleId=obj.sampleId,
                             proportions=obj.proportions
                             )

    @classmethod
    def fromPickle(cls,path) -> ImmuneNetwork:
        with open(path, "rb") as f:
            net_raw = pickle.load(f)
        return cls.fromCustomObject(net_raw)
