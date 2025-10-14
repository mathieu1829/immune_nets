from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.orm import relationship
from sqlalchemy import ForeignKey
from sqlalchemy.types import Date, Float, String, Integer
from sqlalchemy.dialects.postgresql import UUID
import uuid
from datetime import date, datetime
from typing import List
import pandas as pd

from .base import Base
from .networkData import NetworkData

from src.creation.globalSettings import globalSettings
from src.creation.immuneNetwork import ImmuneNetwork

class Network(Base):
    __tablename__ = "network"

    network_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("repertoire.repertoire_id"))
    sample_size: Mapped[int] = mapped_column(Integer, nullable=False)

    name: Mapped[str] =  mapped_column(String(30), nullable = True, default=None)
    algorithm: Mapped[str] =  mapped_column(String(30), nullable = False)
    distance_function: Mapped[str] =  mapped_column(String(30), nullable = False)
    network_algorithm_parameters: Mapped[str] =  mapped_column(String(30), nullable = False)
    version: Mapped[str] =  mapped_column(String(10), nullable = False, default = globalSettings().version)
    username: Mapped[str] =  mapped_column(String(30), nullable = False, default = globalSettings().defaultUsername)
    modificationDate: Mapped[date] =  mapped_column(Date, nullable = False, default=datetime.now())
    creationDate: Mapped[date] =  mapped_column(Date, nullable = False, default=datetime.now())

    source_repertoire: Mapped["Repertoire"] = relationship(back_populates="repertoire_networks") # type: ignore
    network_stats: Mapped[List["NetworkStat"]] = relationship(back_populates="source_network") # type: ignore
    network_edges: Mapped[List["NetworkData"]] = relationship(back_populates="source_network") # type: ignore

    @property
    def algorithmParams(self):
        return eval(self.network_algorithm_parameters)

    def setGraph(self, new_graph):
        self.network_edges = [ 
                           NetworkData(network_id=self.network_id,
                                         r1 = row['r1'],
                                         r2 = row['r2'],
                                        ) 
                           for index,row in new_graph.iterrows()]

    def toImmuneNetwork(self):
        new_graph = pd.DataFrame([{
            "r1": n.r1,
            "r2": n.r2
        } for n in self.network_edges])

        return ImmuneNetwork(graph=new_graph,
                             method=self.algorithm,
                             sampleId=self.repertoire_id,
                             distanceFun=self.distance_function,
                             threshold=eval(self.network_algorithm_parameters)["threshold"],
                             sampleSize=len(self.source_repertoire.clonotypes),
                             name=self.name
                            )

    @classmethod
    def fromImmuneNetwork(cls, network: ImmuneNetwork):
        parameters =  {"threshold":network.threshold}
        parameters = str(parameters)
        new_network = cls(repertoire_id=network.sampleId,
                          algorithm=network.method,
                          distance_function=network.distanceFun,
                          network_algorithm_parameters=parameters,
                          sample_size=network.sampleSize,
                          name = network.name
                          )
        new_network.setGraph(network.graph)
        return new_network
        

