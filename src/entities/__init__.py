from typing import TYPE_CHECKING


__all__ = [
            "ImmuneNetwork",
            "ImmuneRepertoire",
"GraphStats",
            "RepertoireStats",
        ]

if TYPE_CHECKING:
    from .immuneNetwork import ImmuneNetwork
    from .immuneRepertoire import ImmuneRepertoire
    from .graphStats import GraphStats
    from .repertoireStats import RepertoireStats

def __getattr__(name):
    if name == "ImmuneNetwork":
        from .immuneNetwork import ImmuneNetwork
        return ImmuneNetwork
    if name == "ImmuneRepertoire":
        from .immuneRepertoire import ImmuneRepertoire
        return ImmuneRepertoire
    if name == "GraphStats":
        from .graphStats import GraphStats
        return GraphStats
    if name == "RepertoireStats":
        from .repertoireStats import RepertoireStats
        return RepertoireStats

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
