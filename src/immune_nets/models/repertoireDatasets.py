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

class RepertoireDatasets(Base):
    __tablename__ = "REPERTOIRE_DATASETS"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                          primary_key=True,
                                          default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("REPERTOIRES.repertoire_id"), nullable = False)
    dataset_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("DATASETS.dataset_id"), nullable = False)

