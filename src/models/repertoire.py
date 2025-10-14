from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.orm import relationship
from sqlalchemy import ForeignKey
from sqlalchemy.types import Date, Float, String, Integer
from sqlalchemy.dialects.postgresql import UUID
import uuid
from datetime import date, datetime
from typing import List
from src.creation.immuneRepertoire import ImmuneRepertoire
import pandas as pd

from .base import Base
from .clonotypeData import ClonotypeData

from src.creation.globalSettings import globalSettings
from src.creation.algorithms.common_methods import split_tcr_column

class Repertoire(Base):
    __tablename__ = "repertoire"

    repertoire_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                     primary_key=True,
                                                     default=uuid.uuid4)
    name: Mapped[str] =  mapped_column(String(30), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = True)
    version: Mapped[str] =  mapped_column(String(10), nullable = False, default=globalSettings().version)
    username: Mapped[str] =  mapped_column(String(30), nullable = False, default=globalSettings().defaultUsername)
    modificationDate: Mapped[date] =  mapped_column(Date, nullable = False, default=datetime.now())
    creationDate: Mapped[date] =  mapped_column(Date, nullable = False, default=datetime.now())

    datasets: Mapped[List["Dataset"]] = relationship( # type: ignore
        secondary="repertoire_datasets",   
        back_populates="repertoires",
    )
    clonotypes: Mapped[List["ClonotypeData"]] = relationship(back_populates="source_repertoire", cascade="all, delete-orphan") # type: ignore
    
    repertoire_stats: Mapped[List["RepertoireStat"]] = relationship(back_populates="source_repertoire") # type: ignore
    repertoire_networks: Mapped[List["Network"]] = relationship(back_populates="source_repertoire") # type: ignore

    def setClones(self, new_clones: pd.DataFrame):
        self.clonotypes = [ 
                           ClonotypeData(repertoire_id=self.repertoire_id,
                                         proportion = row['proportion'],
                                         tcra_aa = row['tcra_aa'],
                                         tcrb_aa = row['tcrb_aa'],
                                         cdr3s_nt = row['cdr3s_nt'],
                                         inkt_evidence = row['inkt_evidence'],
                                         mait_evidence = row['mait_evidence']
                                        ) 
                           for index,row in new_clones.iterrows()]


    
    
    
