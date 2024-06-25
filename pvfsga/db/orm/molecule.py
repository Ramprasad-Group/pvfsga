from sqlalchemy import Text, Integer, UniqueConstraint, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

from pvfsga.db import Base


class Molecule(Base):
    """Molecules generated from the GA

    Attributes:
        molecule_id: Id of the molecule
        smiles: Smiles string of the molecule
        scaffold_id: Id of the scaffold in the scaffolds table that was used to create
            the molecule
        rgroup_mapping_id: Id of the rgroup_mappings that link which r group was placed
            in which rgroup location in the scaffold
        scaffold: relationship mapping molecule to the scaffold used to create it
    """

    __tablename__ = "molecules"

    molecule_id: Mapped[int] = mapped_column(
        Integer, primary_key=True, autoincrement=True
    )
    smiles: Mapped[str] = mapped_column(Text, nullable=False)
    scaffold_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("scaffolds.scaffold_id"), nullable=False
    )

    scaffold = relationship("Scaffold", back_populates="molecules")
    rgroup_mappings = relationship("RGroupMapping", back_populates="molecule")


class Scaffold(Base):
    """Scaffolds to build molecules from

    Attributes:
        scaffold_id: id of the scaffold in the scaffolds table
        structure: smiles structure of the scaffold to generate molecules from (e.g.,
            C1CC([*:1])CCC1 where [*:1] would be replaced with an R group to generate
            a molecule)
        molecules: Molecules genearted from this structure
        reaction_procedure_mappings: mapping to reactions procedures this scaffold can
            run through

    """

    __tablename__ = "scaffolds"

    scaffold_id: Mapped[str] = mapped_column(
        Integer, primary_key=True, autoincrement=True
    )
    structure: Mapped[str] = mapped_column(Text, nullable=False)

    molecules = relationship("Molecule", back_populates="scaffold")
    reaction_procedure_mappings = relationship("ReactionScaffoldMapping")


class RGroup(Base):
    __tablename__ = "rgroups"

    rgroup_id: Mapped[int] = mapped_column(
        Integer, primary_key=True, autoincrement=True
    )
    smiles: Mapped[str] = mapped_column(Text)


class RGroupMapping(Base):
    """Maps which R groups are in the indexed locations in the scaffolds that created
    the molecule

    Attributes:
        molecule_id: Id in the molecules table that links to the molecule the
            r groups were placed in
        rgroup_id: Id in the rgroups table that links to the rgroup placed in the
            molecule
        index: R group index the rgroup was placed in to from the scaffold
    """

    __tablename__ = "rgroup_mappings"

    molecule_id = mapped_column(
        Integer, ForeignKey("molecules.molecule_id"), nullable=False, primary_key=True
    )
    rgroup_id = mapped_column(
        Integer, ForeignKey("rgroups.rgroup_id"), nullable=False, primary_key=True
    )
    index = mapped_column(Integer, nullable=False, primary_key=True)

    molecule = relationship("Molecule", back_populates="rgroup_mappings")
    rgroup = relationship("RGroup")

    __table_args__ = (
        UniqueConstraint(
            "molecule_id",
            "rgroup_id",
            "index",
            name="unique_rgroup_mapping",
        ),
    )
