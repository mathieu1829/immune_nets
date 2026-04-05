"""
The distributionInference package contains various methods of indentifying the distribution based on probability mass function. All of identifying classes derive from distributionIdentifier interface.
"""

from .distributionIdentifier import DistributionIdentifier
from .klDivergenceIdentifier import KlDivergenceIndentifier
from .kolomogorovIdentifier import KolomogorovIndentifier
from .naiveIdentifier import NaiveIdentifier
