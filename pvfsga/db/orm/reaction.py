from sqlalchemy import ForeignKey, Integer, Text, UniqueConstraint
from sqlalchemy.orm import relationship, Mapped, mapped_column

from pvfsga.db._base import Base


class Reaction(Base):
    """Individual reaction

    Attributes:
        reaction_id: Primary id in reactions table
        smarts: Smarts reaction string (e.g.,
          "[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1>>[*][C:3](=[O:2])[C:5][C:6][C:1][S:7][*]")
    """

    __tablename__ = "reactions"

    reaction_id: Mapped[int] = mapped_column(
        Integer, primary_key=True, autoincrement=True
    )
    smarts: Mapped[str] = mapped_column(Text, nullable=False)


class ReactionStep(Base):
    """Links reaction to specific step in reaction procedure

    Attributes:
        reaction_id (int):
            Id in reaction.Reaction table
        reaction_procedure_id (int):
            Id in reaction.ReactionProcedure table
        step (int):
            Step in the process (indexing starts at 1)
    """

    __tablename__ = "reaction_steps"

    reaction_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("reactions.reaction_id"),
        nullable=False,
        primary_key=True,
    )
    reaction_procedure_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("reaction_procedures.reaction_procedure_id"),
        nullable=False,
        primary_key=True,
    )
    step: Mapped[int] = mapped_column(Integer, nullable=False, primary_key=True)

    reaction = relationship("Reaction")
    procedure = relationship("ReactionProcedure", back_populates="steps")

    __table_args__ = (
        UniqueConstraint(
            "reaction_procedure_id",
            "reaction_id",
            "step",
            name="unique_reaction_step",
        ),
    )


class ReactionProcedure(Base):
    """Polymer reaction that consists of several steps

    Attributes:
        reaction_procedure_id: Id for procedure
        name: name of the reaction procedure
        steps: steps in the procedure
        polymers_mappings: mappings to polymers generated from this procedure
        scaffold_mappings: mappings to scaffolds that can use this procedure
    """

    __tablename__ = "reaction_procedures"

    reaction_procedure_id = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(Text, nullable=False)
    steps = relationship("ReactionStep")
    polymers_mappings = relationship("ReactionPolymerMapping")
    scaffold_mappings = relationship("ReactionScaffoldMapping")


class ReactionScaffoldMapping(Base):
    """Maps reaction procedures to scaffolds that can run through them

    Attributes:
        reaction_procedure_id: id in the reaction_procedures table
        scaffold_id: id in the scaffolds table
    """

    __tablename__ = "reaction_scaffold_mappings"

    reaction_procedure_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("reaction_procedures.reaction_procedure_id"),
        nullable=False,
        primary_key=True,
    )

    scaffold_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("scaffolds.scaffold_id"),
        nullable=False,
        primary_key=True,
    )

    __table_args__ = (
        UniqueConstraint(
            "reaction_procedure_id",
            "scaffold_id",
            name="unique_reaction_scaffold_mapping",
        ),
    )


class ReactionPolymerMapping(Base):
    """Maps reaction procedure, starting molecule, and polymer

    Attributes:
        reaction_procedure_id: Id in reaction.ReactionProcedure table
        pol_id: Id in polymer.Polymer table
        molecule_id: Id in molecule.Molecule table
    """

    __tablename__ = "reaction_polymer_mappings"

    reaction_procedure_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("reaction_procedures.reaction_procedure_id"),
        nullable=False,
        primary_key=True,
    )
    molecule_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("molecules.molecule_id"), nullable=False, primary_key=True
    )
    pol_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("polymers.pol_id"), nullable=False, primary_key=True
    )

    __table_args__ = (
        UniqueConstraint(
            "reaction_procedure_id",
            "pol_id",
            "molecule_id",
            name="unique_reaction_polymer_mapping",
        ),
    )
