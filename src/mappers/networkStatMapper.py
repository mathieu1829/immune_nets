from src.models import NetworkStat, Network
from src.analysis.methods.graphStats import GraphStats
from src.mappers import GraphStatsMapper

class NetworkStatMapper:
    @staticmethod
    def fromGraphStat(stats: GraphStats, network: Network) -> list[NetworkStat]:
        networkStatList = []
        for statName in GraphStats.vectorStatNames() :
            value = getattr(stats,statName)
            networkStat = NetworkStat(network_id=network.network_id,
                                      stat_name=statName,
                                      value=value,
                                      description=""
                                      )
            networkStatList.append(networkStat)
        return networkStatList



    @staticmethod
    def toGraphStat():
        pass
