from typing import Dict, Optional, List

from pydantic import BaseModel


class Molecule(BaseModel):
    """
    Represents a molecule with its structure scaffold, SMILES, and R-groups.

    Attributes:
        structure_scaffold (str):
             The structure scaffold of the molecule.
        smiles (str):
             The SMILES representation of the molecule.
        r_groups (Dict[str, str]):
             Dictionary of R-groups and their values.
    """

    structure_scaffold: str
    smiles: str
    r_groups: Dict[str, str]
    fingerprint: Optional[Dict[str, float]]


class ReactionProcedure(BaseModel):
    """
    Represents the reaction prodecures

    Attributes:
        name (str):
            The name of the reaction procedure
        steps (List[str]):
             a list including all the the reaction steps (smarts)

    """

    name: str
    steps: List[str]


class Polymer(BaseModel):
    """
    Represents a polymer with its Molecule values, polymer smiles, fingerprints,
    and properties.

    Attributes:
        molecule (Molecule):
             The molecule information of the polymer.
        smiles (str):
            The SMILES representation of the polymer.
        canonical_smiles (str):
            The canonicalized SMILES of the polymer.
        fingerprint (Optional[Dict[str, float]]):
             Optional fingerprint of the polymer.
        properties (Optional[Dict[str, float]]):
            Optional properties of the polymer.
    """

    molecule: Molecule
    reaction_proc: ReactionProcedure
    smiles: str
    canonical_smiles: Optional[str]
    fingerprint: Optional[Dict[str, float]]
    properties: Optional[Dict[str, float]]
    fitness: Optional[Dict[str, float]]
