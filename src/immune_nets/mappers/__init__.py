__all__ = [
    "ImmuneNetworkMapper",
    "RepertoireMapper",
    "NetworkMapper",
    "GraphStatsMapper",
    "NetworkStatMapper",
]

def __getattr__(name):
    if name == "ImmuneNetworkMapper":
        from .immuneNetworkMapper import ImmuneNetworkMapper
        return ImmuneNetworkMapper
    if name == "RepertoireMapper":
        from .repertoireMapper import RepertoireMapper
        return RepertoireMapper
    if name == "NetworkMapper":
        from .networkMapper import NetworkMapper
        return NetworkMapper
    if name == "GraphStatsMapper":
        from .graphStatsMapper import GraphStatsMapper
        return GraphStatsMapper
    if name == "NetworkStatMapper":
        from .networkStatMapper import NetworkStatMapper
        return NetworkStatMapper

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
