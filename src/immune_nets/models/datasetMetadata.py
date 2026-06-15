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

class DatasetMetadata(Base): 
    __tablename__ = "DATASET_METADATA"

    dataset_metadata_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                           primary_key=True,
                                                           default=uuid.uuid4)
    dataset_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("DATASETS.dataset_id"))
    name: Mapped[str] =  mapped_column(String(30), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = False)
    value: Mapped[str] =  mapped_column(String, nullable = False)


    source_dataset: Mapped["Dataset"] = relationship(back_populates="dataset_metadata_list") # type: ignore

