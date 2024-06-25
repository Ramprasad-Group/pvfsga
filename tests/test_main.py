import pytest

from pvfsga import main, store_poly_chem
from pvfsga.models import Molecule, Polymer, ReactionProcedure


def fingerprint(molecule: list[Polymer]):
    for che in molecule:
        che.fingerprint = {"test_fp": 1}
    return molecule


def predict(smiles):
    return {"test_property": 1}


@pytest.fixture
def reaction_procedure():
    return ReactionProcedure(
        name="test1",
        steps=[
            "[C:1]1[O:2][C:3](=[S:4])[C:5][C:6]1>>[*][S:4][C:1][C:6][C:5][C:3](=[O:2])[*]"
        ],
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
        r_groups={"R1": "[*]CO"},
        smiles="S=C1C(CO)CCO1",
    )


@pytest.fixture
def silly_polymer_1(silly_chem_1, reaction_procedure):
    return Polymer(
        molecule=silly_chem_1,
        smiles="*C1(COC)C(CO)C(CC)C(=S)O1*",
        reaction_proc=reaction_procedure,
    )


@pytest.fixture
def silly_polymer_2(silly_chem_2, reaction_procedure):
    return Polymer(
        molecule=silly_chem_2,
        smiles="*S=C1C(CO)CCO1*",
        reaction_proc=reaction_procedure,
    )


@pytest.fixture
def list_molecules(silly_chem_1, silly_chem_2):
    return [silly_chem_1, silly_chem_2]


@pytest.fixture
def list_polymers(silly_polymer_1, silly_polymer_2):
    return [silly_polymer_1, silly_polymer_2]


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
        },
        fp_function=fingerprint,
        initial_population_num=10,
        Session=session,
    )


def test_mate(obj_2):
    n = 2
    list_chem = obj_2.mate(
        ["[*]COC", "[*]OO"],
        [
            "C1C([*:1])CC(=S)O1",
        ],
        4,
    )
    assert len(list_chem) == n
    for index, r in enumerate(list_chem):
        for i in range(index + 1, len(list_chem)):
            assert list_chem[index].smiles != list_chem[i].smiles


def test_reproduce(obj_2, list_polymers):
    children = obj_2.reproduce(list_polymers, 4)
    assert len(children) == 4
    assert (
        children[0].molecule.structure_scaffold == "C1([*:3])C([*:1])C([*:2])C(=S)O1"
        or children[0].molecule.structure_scaffold == "S=C1C([*:1])CCO1"
    )


def test_polymerize(obj_2, silly_chem_1):
    pol = obj_2.polymerize(silly_chem_1, obj_2.reaction_procedure_dict["test1"])
    assert pol is not None
    assert "*" in pol.smiles


# Still need to figure out how to append fingerprints
def test_polymerize_all(list_molecules, obj_2):
    store_poly_chem.store_molecules(obj_2, list_molecules)
    polymer_list = obj_2.polymerize_all(list_molecules)
    for poly in polymer_list:
        assert poly.smiles is not None
        assert poly.fingerprint is None


def test_fix(silly_chem_1):
    silly_chem_1_sub = silly_chem_1
    assert silly_chem_1_sub.structure_scaffold == "C1([*:3])C([*:1])C([*:2])C(=S)O1"
