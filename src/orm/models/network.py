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

class Network(Base):
    __tablename__ = "network"

    network_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("repertoire.repertoire_id"))
    sample_size: Mapped[int] = mapped_column(Integer, nullable=False)

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

    @property
    def graph(self):
        if not hasattr(self, "_clones_df"):
            self._graph_df = pd.DataFrame([{
                "r1": n.r1,
                "r2": n.r2
            } for n in self.network_edges])
        if self._graph_df.empty:
            self._graph_df = pd.DataFrame({"r1":[],"r2":[]})

        return self._graph_df

    @graph.setter
    def graph(self, new_graph):
        self._graph_df = new_graph
        self.network_edges = [ 
                           NetworkData(network_id=self.network_id,
                                         r1 = row['r1'],
                                         r2 = row['r2'],
                                        ) 
                           for index,row in self._graph_df.iterrows()]

    @classmethod
    def createFullNetwork(cls, repertoire_id, graph, method, distanceFun, threshold, sampleSize):
        parameters =  {"threshold":threshold}
        parameters = str(parameters)
        network = Network(repertoire_id=repertoire_id,
                          algorithm=method,
                          distance_function=distanceFun,
                          network_algorithm_parameters=parameters,
                          sample_size=sampleSize
                          )
        network.graph = graph
        return network
        

