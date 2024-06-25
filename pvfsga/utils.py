from itertools import product
import logging

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from pvfsga.models import Molecule

logger = logging.getLogger(__name__)


def create_all_molecule_list(
    r_groups: list[str], scaffolds: list[str]
) -> list[Molecule]:
    """Create a list of all possible molecules itetation through all
     repeating combination of the R-groups

    Checks for any errors in R-group addition and then
    logger.debugs a statement declaring if it found one

    O(n^{r_loc})

    Returns:
        a list containing the Molecules
    """
    molecules = []
    for scaffold in scaffolds:
        r_num_positions = scaffold.count("[*:")
        tot_combinations = len(r_groups) ** r_num_positions
        if tot_combinations > 1000000:
            logger.debug(f"too many molecules {tot_combinations}, will not preceed")
            return []

        perm_list = list(product(r_groups, repeat=r_num_positions))
        for perm in perm_list:
            mol = create_molecule(scaffold, perm)
            if mol is not None:
                molecules.append(mol)

    return molecules


def create_molecule(scaffold: str, r_groups: list[str]):
    """Create molecule from r_group list

    Requires r_group list is same length as number of viable locations in sub
    """
    r_num_positions = scaffold.count("[*:")
    if r_num_positions != len(r_groups):
        raise ValueError(
            f"Number of r_groups in {r_groups} must equal number of sub locations in"
            + f"{scaffold}"
        )
    mol = scaffold
    for i, R in enumerate(r_groups):
        mol = insert_r_group(mol, i + 1, R)
        if mol is None:
            logger.info("The following molecule could not be processed", mol)
            return None
    r_groups_dict = {f"R{i+1}": r_groups[i] for i in range(r_num_positions)}
    return Molecule(
        smiles=mol,
        structure_scaffold=scaffold,
        r_groups=r_groups_dict,
    )


def insert_r_group(scaffold: str, r_pos_index: int, R: str):
    r_string = "[*:" + str(r_pos_index) + "]"
    if R == "[*][H]" or R == "":
        return scaffold.replace(r_string, "[H]")
    mR = Chem.MolFromSmiles(R.replace("[*]", "[Bi]"))
    Chem.AddHs(mR, explicitOnly=True)
    mscaffold = scaffold.replace(r_string, "[Bi]")
    mscaffold = Chem.MolFromSmiles(mscaffold)
    rxn = rdChemReactions.ReactionFromSmarts(
        "[*:1][Bi:2].[*:3][Bi:4]>>[*:1][*:3].[Bi:2].[Bi:4]"
    )
    try:
        smiles = Chem.MolToSmiles(rxn.RunReactants((mscaffold, mR))[0][0])
    except Exception as e:
        logger.warning(f"Could not insert {R} group in {scaffold} due to {e}")
        smiles = None
    return smiles
