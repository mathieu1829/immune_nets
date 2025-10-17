import pandas as pd

from src.analysis.methods.graphStats import GraphStats

class GraphStatsMapper:
    @staticmethod
    def toStatVector(stats:GraphStats):
        return [ 
                float(stats.isolateVerticeRatio),
                float(stats.edgeDensity),
                float(stats.density),
                float(stats.meanEccentricity),
                float(stats.giantComponent),
                float(stats.meanDegree),
                float(stats.componentCount),
                float(stats.meanComponentSize)
               ]

    @staticmethod
    def toList(stats:GraphStats):
        return [ 
                stats.isolateVerticeRatio,
                stats.edgeDensity,
                stats.density,
                stats.eccentricity,
                stats.giantComponent,
                stats.degreeDistribution,
                stats.componentCount,
                stats.componentSizeDistribution
               ]

    @staticmethod
    def colnames():
        return [
                "isolatedVerticeRatio",
                "edgeDensity",
                "density",
                "meanEccentricity",
                "giantComponent",
                "meanDegree",
                "componentCount",
                "meanComponentSize"
               ]

    @classmethod
    def toDataframe(cls, statList: list[GraphStats]) -> pd.DataFrame:
        data = [ cls.toList(stat) for stat in statList ]
        df = pd.DataFrame(data=data, columns=pd.Index(cls.colnames()))
        return df
        

