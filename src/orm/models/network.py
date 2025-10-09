from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.orm import relationship
from sqlalchemy import ForeignKey
from sqlalchemy.types import Date, Float, String, Integer
from sqlalchemy.dialects.postgresql import UUID
import uuid
from datetime import date,time
from typing import List

from .base import Base

class Network(Base):
    __tablename__ = "network"

    network_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("repertoire.repertoire_id"))

    algorithm: Mapped[str] =  mapped_column(String(30), nullable = False)
    distance_function: Mapped[str] =  mapped_column(String(30), nullable = False)
    network_algorithm_parameters: Mapped[str] =  mapped_column(String(30), nullable = False)
    version: Mapped[str] =  mapped_column(String(10), nullable = False)
    username: Mapped[str] =  mapped_column(String(30), nullable = False)
    modificationDate: Mapped[date] =  mapped_column(Date, nullable = False)
    creationDate: Mapped[date] =  mapped_column(Date, nullable = False)

    source_repertoire: Mapped["Repertoire"] = relationship(back_populates="repertoire_networks") # type: ignore
    network_stats: Mapped[List["NetworkStat"]] = relationship(back_populates="source_network") # type: ignore
    network_edges: Mapped[List["NetworkData"]] = relationship(back_populates="source_network") # type: ignore

