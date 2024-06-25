from typing import Optional

from sqlalchemy import Integer, ForeignKey, Text, Float, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pvfsga.db import Base


class Polymer(Base):
    """

    Attributes:
        pol_id: id of the generated polymer
        smiles: polymer smiles string
        canonical_smiles: canonical version of the smiles string
        fingerprint: dictionary of the fingerprint stored in the db
        reaction_procedure_mappings: mapping to table that links to reaction procedures
            and molecules that can be used to create this polymer
    """

    __tablename__ = "polymers"

    pol_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    smiles: Mapped[str] = mapped_column(Text, nullable=False)
    canonical_smiles: Mapped[Optional[str]] = mapped_column(Text)

    reaction_procedure_mappings = relationship("ReactionPolymerMapping")
    properties = relationship("PolymerProperty")
    fingerprint = relationship("PolymerFingerprint")


class PolymerProperty(Base):
    __tablename__ = "polymer_properties"

    pol_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("polymers.pol_id"), nullable=False, primary_key=True
    )
    prop: Mapped[str] = mapped_column(Text, nullable=False, primary_key=True)
    prop_unit: Mapped[str] = mapped_column(Text)
    value: Mapped[int] = mapped_column(Float, nullable=False)

    __table_args__ = (
        UniqueConstraint(
            "pol_id",
            "prop",
            name="unique_polymer_property",
        ),
    )


class PolymerFingerprint(Base):
    __tablename__ = "polymer_fingerprints"

    pol_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("polymers.pol_id"), nullable=False, primary_key=True
    )
    key: Mapped[str] = mapped_column(Text, nullable=False, primary_key=True)
    value: Mapped[int] = mapped_column(Float, nullable=False)

    __table_args__ = (
        UniqueConstraint(
            "pol_id",
            "key",
            name="unique_polymer_fingerprint",
        ),
    )


class PolymerGeneration(Base):
    """There can technically be the same polymer multiple times in a generation"""

    __tablename__ = "polymer_generation"

    # needed because multiple of the same polymer can be made in a generation
    gen_id: Mapped[int] = mapped_column(Integer, nullable=False, primary_key=True)
    pol_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("polymers.pol_id"), nullable=False
    )
    generation: Mapped[int] = mapped_column(Integer, nullable=False)
