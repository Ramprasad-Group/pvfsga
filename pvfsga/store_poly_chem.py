from typing import List
import sys

from pvfsga.db import adder
from pvfsga.models import Molecule, Polymer


def store_molecules(obj_GA, molecule_list_to_store: List[Molecule]):
    """
    Stores molecules' information in the database.

    Args:
        molecule_list_to_store (list[Molecule]): A list of Molecule objects to store.

    Returns:
        None
    """
    with obj_GA.Session() as session:
        for final_chem in molecule_list_to_store:
            if final_chem.smiles not in obj_GA.molecules_dict:
                chem_id = adder.add_molecule_and_get_id(
                    final_chem.smiles,
                    obj_GA.molecule_scaffold_dict[final_chem.structure_scaffold],
                    session,
                )

                obj_GA.molecules_dict[final_chem.smiles] = chem_id
                r_group_ind = 1
                for key, value in final_chem.r_groups.items():
                    adder.add_rgroup_mapping(
                        chem_id, obj_GA.r_group_dict[value], r_group_ind, session
                    )
                    r_group_ind = r_group_ind + 1
        session.commit()


def store_polymer(obj_GA, polymer_list_to_store: List[Polymer]):
    """
    Stores polymer information in the database.

    Args:
        polymer_list_to_store (list[Polymer]): A list of Polymer objects to store.

    Returns:
        None
    """
    obj_GA.polymer_dict_current_gen = {}
    with obj_GA.Session() as session:
        prog = 0
        for final_chem in polymer_list_to_store:
            if obj_GA.verbose > 0:
                PROG = 100 * round(prog / (len(polymer_list_to_store)), 3)
                sys.stdout.write("storing \r%d%% done" % PROG)
                sys.stdout.flush()
                prog += 1
            if (
                final_chem.canonical_smiles not in obj_GA.polymer_dict
                and final_chem.fingerprint is not None
            ):
                pol_id = adder.add_polymer_and_get_id(
                    final_chem.smiles, final_chem.canonical_smiles, session
                )
                obj_GA.polymer_dict_current_gen[final_chem.canonical_smiles] = pol_id
                obj_GA.polymer_dict[final_chem.canonical_smiles] = pol_id
                obj_GA.polymer_fingerprint[
                    final_chem.canonical_smiles
                ] = final_chem.fingerprint
                for key, value in final_chem.fingerprint.items():
                    adder.add_polymer_fingerprint(key, value, pol_id, session)

                adder.add_reaction_polymer_mapping(
                    obj_GA.procedure_dict[final_chem.reaction_proc.name],
                    obj_GA.molecules_dict[final_chem.molecule.smiles],
                    pol_id,
                    session,
                )
            adder.add_polymer_generation(
                obj_GA.generation,
                obj_GA.polymer_dict[final_chem.canonical_smiles],
                session,
            )
            obj_GA.polymer_dict_current_gen[
                final_chem.canonical_smiles
            ] = obj_GA.polymer_dict[final_chem.canonical_smiles]

        session.commit()
