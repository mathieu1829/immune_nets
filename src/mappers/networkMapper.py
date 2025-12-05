import pandas as pd
import numpy as np
import uuid

from src.entities import ImmuneNetwork
from src.models import Network

class NetworkMapper:
    @staticmethod
    def toImmuneNetwork(network: Network) -> ImmuneNetwork:
        new_graph = pd.DataFrame([{
            "r1": n.r1,
            "r2": n.r2
        } for n in network.network_edges])

        return ImmuneNetwork(graph=new_graph,
                             method=network.algorithm,
                             sampleId=network.repertoire_id,
                             distanceFun=network.distance_function,
                             threshold=eval(network.network_algorithm_parameters)["threshold"],
                             sampleSize=len(network.source_repertoire.clonotypes),
                             proportions=np.array([ clone.proportion  for clone in network.source_repertoire.clonotypes]),
                             name=network.name
                            )

    @staticmethod
    def fromImmuneNetwork( network: ImmuneNetwork) -> Network:
        parameters =  {"threshold":network.threshold}
        parameters = str(parameters)
        new_network = Network(network_id=network.network_id,
                          repertoire_id=network.sampleId,
                          algorithm=network.method,
                          distance_function=network.distanceFun,
                          network_algorithm_parameters=parameters,
                          name = network.name
                          )
        new_network.setGraph(network.graph)
        return new_network
