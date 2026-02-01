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

from src.globalSettings import GlobalSettings
from src.entities import ImmuneNetwork

class Network(Base):
    __tablename__ = "network"

    network_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("repertoire.repertoire_id"))

    name: Mapped[str] =  mapped_column(String(30), nullable = True, default=None)
    algorithm: Mapped[str] =  mapped_column(String(50), nullable = False)
    distance_function: Mapped[str] =  mapped_column(String(50), nullable = False)
    network_algorithm_parameters: Mapped[str] =  mapped_column(String(50), nullable = False)
    version: Mapped[str] =  mapped_column(String(10), nullable = False, default = GlobalSettings().version)
    username: Mapped[str] =  mapped_column(String(30), nullable = False, default = GlobalSettings().defaultUsername)
    modificationDate: Mapped[date] =  mapped_column(Date, nullable = False, default=datetime.now())
    creationDate: Mapped[date] =  mapped_column(Date, nullable = False, default=datetime.now())

    source_repertoire: Mapped["Repertoire"] = relationship(back_populates="repertoire_networks") # type: ignore
    network_stats: Mapped[List["NetworkStat"]] = relationship(back_populates="source_network", cascade="all, delete-orphan") # type: ignore
    network_edges: Mapped[List["NetworkData"]] = relationship(back_populates="source_network", cascade="all, delete-orphan") # type: ignore

    @property
    def algorithmParams(self):
        return eval(self.network_algorithm_parameters)

    def setGraph(self, new_graph: pd.DataFrame):
        new_graph = new_graph.dropna(subset=["r1", "r2"])

        new_graph["r1"] = new_graph["r1"].astype(int)
        new_graph["r2"] = new_graph["r2"].astype(int)

        self.network_edges = [ 
                           NetworkData(network_id=self.network_id,
                                         r1 = int(row['r1']),
                                         r2 = int(row['r2']),
                                        ) 
                           for index,row in new_graph.iterrows()]



        

