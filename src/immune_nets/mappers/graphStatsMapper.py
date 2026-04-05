import pandas as pd

from immune_nets.entities import GraphStats

class GraphStatsMapper:
    @classmethod
    def toDataframe(cls, statList: list[GraphStats]) -> pd.DataFrame:
        data = [ stat.toList() for stat in statList ]
        df = pd.DataFrame(data=data, columns=pd.Index(GraphStats.listStatNames()))
        return df
        

