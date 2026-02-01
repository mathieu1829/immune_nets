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

class NetworkData(Base):
    __tablename__ = "network_data"

    network_record_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    network_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("network.network_id"), nullable = False)
    r1: Mapped[int] = mapped_column(Integer, nullable = False)
    r2: Mapped[int] = mapped_column(Integer, nullable = False)
    weights: Mapped[str] = mapped_column(String, nullable = True, default = None)
    r1_weights: Mapped[str] = mapped_column(String, nullable = True, default = None)
    r2_weights: Mapped[str] = mapped_column(String, nullable = True, default = None)


    source_network: Mapped["Network"] = relationship(back_populates="network_edges") # type: ignore

