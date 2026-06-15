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

class ClonotypeData(Base): 
    __tablename__ = "CLONOTYPE_DATA"

    clone_record_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                            primary_key=True,
                                                            default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("REPERTOIRES.repertoire_id"), nullable=False)
    proportion: Mapped[float] = mapped_column(Float, nullable = False)
    tcra_aa: Mapped[str] = mapped_column(String, nullable = True)
    tcrb_aa: Mapped[str] = mapped_column(String, nullable = True)
    cdr3s_nt: Mapped[str] = mapped_column(String, nullable = True)
    inkt_evidence: Mapped[str] = mapped_column(String, nullable = True)
    mait_evidence: Mapped[str] = mapped_column(String, nullable = True)

    source_repertoire: Mapped[List["Repertoire"]] = relationship(back_populates="clonotypes") # type: ignore
