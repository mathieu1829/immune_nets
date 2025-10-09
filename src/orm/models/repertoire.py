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

class Repertoire(Base):
    __tablename__ = "repertoire"

    repertoire_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                     primary_key=True,
                                                     default=uuid.uuid4)
    name: Mapped[str] =  mapped_column(String(30), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = False)
    version: Mapped[str] =  mapped_column(String(10), nullable = False)
    username: Mapped[str] =  mapped_column(String(30), nullable = False)
    modificationDate: Mapped[date] =  mapped_column(Date, nullable = False)
    creationDate: Mapped[date] =  mapped_column(Date, nullable = False)

    datasets: Mapped[List["Dataset"]] = relationship( # type: ignore
        secondary="repertoire_datasets",   
        back_populates="repertoires",
    )
    clonotypes: Mapped[List["ClonotypeData"]] = relationship(back_populates="source_repertoire") # type: ignore
    
    repertoire_stats: Mapped[List["RepertoireStat"]] = relationship(back_populates="source_repertoire") # type: ignore
    repertoire_networks: Mapped[List["Network"]] = relationship(back_populates="source_repertoire") # type: ignore

