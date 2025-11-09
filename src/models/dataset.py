from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.orm import relationship
from sqlalchemy import ForeignKey
from sqlalchemy.types import Date, Float, String, Integer
from sqlalchemy.dialects.postgresql import UUID
import uuid
from datetime import date,datetime
from typing import List

from .base import Base

from src.globalSettings import GlobalSettings

class Dataset(Base):
    __tablename__ = "dataset"

    dataset_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                  primary_key=True,
                                                  default=uuid.uuid4)
    name: Mapped[str] =  mapped_column(String(30), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = False)
    version: Mapped[str] =  mapped_column(String(10), nullable = False, default=GlobalSettings().version)
    username: Mapped[str] =  mapped_column(String(30), nullable = False, default=GlobalSettings().defaultUsername)
    modificationDate: Mapped[date] =  mapped_column(Date, nullable = False, default=datetime.now())
    creationDate: Mapped[date] =  mapped_column(Date, nullable = False, default=datetime.now())

    dataset_metadata_list: Mapped[List["DatasetMetadata"]] = relationship(back_populates="source_dataset", cascade="all, delete-orphan") # type: ignore
    repertoires: Mapped[List["Repertoire"]] = relationship( # type: ignore
        secondary="repertoire_datasets",   
        back_populates="datasets",
    )
