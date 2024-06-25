from typing import List
import pytest
from canonicalize_psmiles.canonicalize import canonicalize
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from pvfsga import main
from pvfsga.models import Molecule, Polymer, ReactionProcedure

r_groups = [
    "[*][H]",
    "[*]CO",
    "[*]O",
    "[*]COO",
    "[*]COOO",
    "[*]COC",
    "[*]C=C",
    "[*]CC",
    "[*][CH3]",
    "[*]OO",
    "[*]C",
    "[*]N=O",
    "[*]OCl",
    "[*]SO",
]


def fp(polymers: List[Polymer]) -> List[Polymer]:
    for polymer in polymers:
        polymer.fingerprint = {"fp1": 1}
    return polymers


@pytest.fixture
def reaction_procedure():
    return ReactionProcedure(
        name="test1",
        steps=[
            "[C:1]1[O:2][C:3](=[S:4])[C:5][C:6]1>>[*][S:4][C:1][C:6][C:5][C:3](=[O:2])[*]"
        ],
    )


@pytest.fixture
def obj_2(session):
    return main.GA_society(
        ["[*]COC", "[*]OO", "[*]CO", "[*]CC"],
        {
            "test1": [
                "[C:1]1[O:2][C:3](=[S:4])[C:5][C:6]1>>[*][S:4][C:1][C:6][C:5][C:3](=[O:2])[*]"
            ]
        },
        {
            "C1C([*:1])CC(=S)O1": ["test1"],
            "C1([*:3])C([*:1])C([*:2])C(=S)O1": ["test1"],
            "S=C1C([*:1])CCO1": ["test1"],
            "CCCC1C([*:1])CC(=S)O1": ["test1"],
        },
        fp_function=fp,
        initial_population_num=10,
        Session=session,
    )


@pytest.fixture
def silly_chem_1():
    return Molecule(
        structure_scaffold="C1([*:3])C([*:1])C([*:2])C(=S)O1",
        r_groups={
            "R1": "[*]CO",
            "R2": "[*]CC",
            "R3": "[*]COC",
        },
        smiles="C1(COC)C(CO)C(CC)C(=S)O1",
    )


@pytest.fixture
def silly_chem_2():
    return Molecule(
        structure_scaffold="S=C1C([*:1])CCO1",
        r_groups={"R1": "[*]OO"},
        smiles="S=C1C(OO)CCO1",
    )


@pytest.fixture
def silly_chem_3():
    return Molecule(
        structure_scaffold="S=C1C([*:1])CCO1",
        r_groups={"R1": "[*]COC"},
        smiles="S=C1C(COC)CCO1",
    )


@pytest.fixture
def silly_polymer_1(silly_chem_1, reaction_procedure):
    rxn = rdChemReactions.ReactionFromSmarts(reaction_procedure.steps[0])
    polymer = Chem.MolToSmiles(
        rxn.RunReactants((Chem.MolFromSmiles(silly_chem_1.smiles),))[0][0]
    )
    return Polymer(
        molecule=silly_chem_1,
        smiles=polymer,
        canonical_smiles=canonicalize(polymer),
        reaction_proc=reaction_procedure,
    )


@pytest.fixture
def silly_polymer_2(silly_chem_2, reaction_procedure):
    rxn = rdChemReactions.ReactionFromSmarts(reaction_procedure.steps[0])
    polymer = Chem.MolToSmiles(
        rxn.RunReactants((Chem.MolFromSmiles(silly_chem_2.smiles),))[0][0]
    )
    return Polymer(
        molecule=silly_chem_2,
        smiles=polymer,
        canonical_smiles=canonicalize(polymer),
        reaction_proc=reaction_procedure,
    )


@pytest.fixture
def silly_polymer_3(silly_chem_3, reaction_procedure):
    rxn = rdChemReactions.ReactionFromSmarts(reaction_procedure.steps[0])
    polymer = Chem.MolToSmiles(
        rxn.RunReactants((Chem.MolFromSmiles(silly_chem_3.smiles),))[0][0]
    )
    return Polymer(
        molecule=silly_chem_3,
        smiles=polymer,
        canonical_smiles=canonicalize(polymer),
        reaction_proc=reaction_procedure,
    )


@pytest.fixture
def list_molecules(silly_chem_1, silly_chem_2):
    return [silly_chem_1, silly_chem_2]


@pytest.fixture
def list_polymers(silly_polymer_1, silly_polymer_2, silly_polymer_3):
    return [silly_polymer_1, silly_polymer_2, silly_polymer_3]


# Error on purpose, to see see whether kids are allowed to be exactly like their
# parents or not
def test_crossover(obj_2, silly_chem_1, silly_chem_2):
    list_mate = obj_2.crossover(silly_chem_1, silly_chem_2, 4)
    assert len(list_mate) == 4


def test_crossover_all(list_polymers, obj_2):
    total_list_mate = obj_2.crossover_all(list_polymers)
    assert len(total_list_mate) == 8
