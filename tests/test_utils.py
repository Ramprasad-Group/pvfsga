from pvfsga.utils import create_all_molecule_list


def test_create_all_molecule_list():
    r_groups = ["[*]C", "[*]O", "[*]F"]
    scaffolds = ["C1([*:3])C([*:1])C([*:2])C(=S)O1"]
    assert len(create_all_molecule_list(r_groups, scaffolds)) == 3**3
