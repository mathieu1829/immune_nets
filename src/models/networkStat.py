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

class NetworkStat(Base):
    __tablename__ = "network_stat"

    network_stat_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    network_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("network.network_id"), nullable = False)
    stat_name: Mapped[str] =  mapped_column(String(30), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = False)
    value: Mapped[str] =  mapped_column(String(120), nullable = False)

    source_network: Mapped["Network"] = relationship(back_populates="network_stats") # type: ignore

