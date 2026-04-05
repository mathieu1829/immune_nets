import re
import numpy as np
import ast
import igraph as ig

from sqlalchemy.orm import Session

from immune_nets.models import NetworkStat, Network
from immune_nets.entities import ImmuneNetwork
from immune_nets.entities import GraphStats
from immune_nets.factories import ImmuneNetworkFactory
from immune_nets.mappers import NetworkMapper
from immune_nets.db import engine

class NetworkStatMapper:
    @classmethod
    def fromGraphStat(cls, stats: GraphStats, network: Network | ImmuneNetwork) -> list[NetworkStat]:
        networkStatList = []
        for statName, statValue in stats.__dict__.items() :
            if statName == "graph":
                continue
            if statName == "components":
                continue

            dbDatatype = cls.inferType(statName, statValue)
            match dbDatatype:
                case "ndarray":
                    dbValue = str(statValue.tolist())
                case _:
                    dbValue = str(statValue)

                

            networkStat = NetworkStat(network_id=network.network_id,
                                      stat_name=statName,
                                      value=dbValue,
                                      datatype=dbDatatype
                                      )
            networkStatList.append(networkStat)
        return networkStatList



    @staticmethod
    def toGraphStat(networkStatList) -> GraphStats:
        baseStats = GraphStats(ImmuneNetworkFactory.createEmpty())
        for networkStat in networkStatList:
            if networkStat.datatype == "ndarray":
                processedValue = np.array(ast.literal_eval(networkStat.value))
            else:
                processedValue = ast.literal_eval(networkStat.value)
            
            setattr(baseStats,networkStat.stat_name,processedValue)
        with Session(engine) as session:
            sourceNetwork: Network|None = session.get(Network, networkStatList[0].network_id)
            if sourceNetwork:
                immuneNet: ImmuneNetwork = NetworkMapper.toImmuneNetwork(sourceNetwork)
            else:
                raise ValueError(f"No network with id: {networkStatList[0].network_id}")
            baseStats.graph = ig.Graph(immuneNet.graph.to_numpy())
            baseStats.components = baseStats.graph.components()
        return baseStats

    @staticmethod
    def inferType(name: str, value):
        if re.search(r".*Distribution",name):
            return "Dict[str,float]"
        elif re.search(r".*List",name):
            return f"List[{type(value[0]).__name__}]"
        else:
            return type(value).__name__

