import pytest

from pvfsga import models


@pytest.fixture
def molecule():
    return models.Molecule(
        structure_scaffold="CCCC1C([*:1])CC(=S)O1",
        smiles="CCCC1C(CO)CC(=S)O1",
        r_groups={"R1": "CO"},
    )


@pytest.fixture
def reaction_procedure():
    return models.ReactionProcedure(
        name="test1",
        steps=[
            "[C:1]1[O:2][C:3](=[S:4])[C:5][C:6]1>>[*][S:4][C:1][C:6][C:5][C:3](=[O:2])[*]"
        ],
    )


def test_molecule():
    chem = models.Molecule(
        structure_scaffold="CCCC1C([*:1])CC(=S)O1",
        smiles="CCCC1C(CO)CC(=S)O1",
        r_groups={"R1": "CO"},
    )
    assert chem.smiles == "CCCC1C(CO)CC(=S)O1"
    assert chem.structure_scaffold == "CCCC1C([*:1])CC(=S)O1"
    assert "R1" in chem.r_groups


def test_polymer(molecule, reaction_procedure):
    pol = models.Polymer(
        molecule=molecule,
        smiles="test",
        fingerprint={"test": 1},
        properties=None,
        reaction_proc=reaction_procedure,
    )
    assert pol.smiles == "test"
    assert pol.molecule.smiles == "CCCC1C(CO)CC(=S)O1"
    pol.molecule.smiles = "test"
    assert molecule.smiles == "test"
