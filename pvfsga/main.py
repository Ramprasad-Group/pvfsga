from typing import Optional
import logging
from datetime import datetime
import sys

from rdkit import Chem
from rdkit.Chem import rdChemReactions
import numpy as np
from sqlalchemy.orm.session import Session
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from canonicalize_psmiles.canonicalize import canonicalize

from pvfsga.db import Base
from pvfsga.db import adder
from pvfsga.models import Molecule, Polymer, ReactionProcedure
from pvfsga import store_poly_chem
from pvfsga import utils


logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.propagate = False


class GA_society:
    """GA_society is the class that represents the population of molecules
        and the generative algorithm attributes/functions.

    Attributes:
        r_groups (list[str]):
            the R-groups that can be used for the scaffolds

        scaffolds (list[str]):
            a list of all base scaffolds used to build new molecules

        fp_function (callable):
            Fingerprinting fucntion that takes in a list[Polymer], and return
            list[Polymer] with Polymer.fingerprint

        initial_population ( list[Molecule]):
             list of Molecule containing initial population. If
            initial_population is None , random initial population is used.
            Defualt None

        initial_population_num (int):
            number of species in the initial population

        current_population (list):
            A list of the current population, that corresponds to a
            certain generation number

        num_families (int):
            Number of families that will propagate. Default 50

        num_children_per_family (int):
            Number of children per pair of parents. Default is 3.

        generation (int):
            Current generation of polymers/molecules.

        partner_selection (str):
            str representing how parents choose their mate.
            'diversity' means highest scoring parents choose
            partner based on least similar tanimoto similarity
            score. 'random' means partner chosen randomly.
            Default 'diversity'. 'all' means all parents are mated together

        random_seed (int):
            Random seed to use for nation. If None, no
            random seed is used.

        rng (numpy.random.Generator):
            The random number generator.

        polymer_dict (dict[str:int]):
            dictionary mapping each polymer to its polymer_id of the db

        molecules_dict (dict[str:int]):
            dictionary mapping each molecule to its chem_id of the db

        r_group_dict (dict[str:int]):
            dictionary mapping each r_group to its r_group_id of the db

        reaction_steps_dict (dict[str:int]):
            dictionary mapping each reaction step to its reach step id of the db

        molecule_scaffold_dict (dict[str:int]):
            dictionary mapping each chemica; scaffold to its id at the db

        procedure_dict(dict[str:int]):
            dictionary mapping each reaction procedure to its id at the db

        scaffold_reaction_dict(dict[str:list[str]]):
            maps each scaffold/ to a list of all the viable reaction procedures

        polymer_fingerprint(dict[str:dict]):
            maps each polymer to its fingerprint dictionary

        polymer_property(dict[str:dict]):
            maps each polymer to its property dictionary


        Session(Optional[Session]):
            DB session connection

    """

    def __init__(
        self,
        r_groups: list[str],
        reaction_procedure_dict: dict[str, list[str]],
        scaffold_reaction_dict: dict[str, list[str]],
        fp_function: callable,
        initial_population: list[Molecule] = None,
        create_random: bool = True,
        initial_population_num: int = 20,
        num_families: int = 50,
        num_children_per_family: int = 3,
        max_molecule_atom_count: int = 100,
        generation: int = 1,
        partner_selection: str = "all",
        random_seed: Optional[int] = None,
        Session: Optional[Session] = None,
        verbose: int = 0,
    ):
        """
        Args:
            reaction_procedure_dict(dict[str,list[str]]):
                Key holds the name of the reaction_procedures
                values holds a list of the reaction steps (smarts)

        """
        self.verbose = verbose
        self.r_groups = r_groups
        self.scaffold_reaction_dict = scaffold_reaction_dict
        self.scaffolds = list(scaffold_reaction_dict.keys())
        self.polymer_dict = {}
        self.polymer_fingerprint = {}
        self.molecules_dict = {}
        self.r_group_dict = {}
        self.reaction_steps_dict = {}
        self.molecule_scaffold_dict = {}
        self.procedure_dict = {}
        self.polymer_property = {}
        self.max_molecule_atom_count = max_molecule_atom_count

        self.scaffolds = list(set(self.scaffolds))
        if len(self.scaffolds) != len(list(scaffold_reaction_dict.keys())):
            logger.info("A duplicte  was passed. Removed from the set")
        for sub in self.scaffolds:
            r_num_positions = sub.count("[*:")
            if r_num_positions == 0:
                self.scaffolds.remove(sub)
                logger.info(
                    "A  with no positions for R-groups was passed. Removed from the set"
                )

        self.reaction_procedure_dict = {}
        for key in reaction_procedure_dict:
            self.reaction_procedure_dict[key] = ReactionProcedure(
                name=key, steps=reaction_procedure_dict[key]
            )

        if Session is None:
            db_server = (
                f"sqlite:///run_{datetime.now().strftime('%Y_%m_%d_%H-%M-%S')}.sqlite"
            )
            engine = create_engine(db_server)
            Base.metadata.create_all(engine)
            self.Session = sessionmaker(engine)
        else:
            self.Session = Session

        with self.Session() as session:
            for procedure in self.reaction_procedure_dict.values():
                procedure_id = adder.add_reaction_procedure_and_or_get_id(
                    procedure, session
                )
                self.procedure_dict[procedure.name] = procedure_id

            for rgroup in self.r_groups:
                rgroup_id = adder.add_rgroup_and_or_get_id(rgroup, session)
                self.r_group_dict[rgroup] = rgroup_id

            for scaffold in self.scaffolds:
                scaffold_id = adder.add_scaffold_and_or_get_id(scaffold, session)

                self.molecule_scaffold_dict[scaffold] = scaffold_id
                for reaction_proc_name in self.scaffold_reaction_dict[scaffold]:
                    adder.add_reaction_scaffold_mappings(
                        scaffold_id, self.procedure_dict[reaction_proc_name], session
                    )

            session.commit()

        self.fp_function = fp_function
        self.rng = np.random.default_rng(seed=random_seed)

        if initial_population is None:
            if create_random:
                self.initial_population = self.create_random_intital_molecule_list(
                    self.r_groups,
                    self.scaffolds,
                    initial_population_num,
                )
            else:
                 self.initial_population  = utils.create_all_molecule_list(                    
                    self.r_groups,
                    self.scaffolds
                )
        else:
            self.initial_population = initial_population

        self.current_population = self.initial_population
        store_poly_chem.store_molecules(
            self, [x.copy() for x in self.initial_population]
        )

        self.num_families = num_families
        self.num_children_per_family = num_children_per_family
        self.generation = generation
        self.partner_selection = partner_selection

    def create_random_intital_molecule_list(
        self, r_groups: list[str], scaffolds: list[str], num_mols: int
    ) -> list[Molecule]:
        """Creates initial population pydantic containing num_mols molecule
        chosen randomly from the r-groups and scaffolds

        Checks for any errors in R-group addition and then
        logger.debugs a statement declaring if it found one

        Args:
            r_groups (list[str]):
                list of all R-groups
            scaffolds (list[str]):
                list of all scaffolds
            num_mols (int):
                The number of Molecules to create

        Returns:
            chem_list(list[Molecule]):
                a list containing the Molecule values.
        """

        molecules = []  # stores all the newly generated molecules
        generated_molecules = set()  # ensure we only create a new molecule once
        scaffold_id, failed_attempts = 0, 0
        while len(molecules) < num_mols and failed_attempts < 1000:
            scaffold = scaffolds[scaffold_id]
            prog = 100 * round(len(molecules) / (num_mols), 3)
            if self.verbose > 0:
                print(f"initial population {prog} % done", end="\r")

            r_num_positions = scaffold.count("*:")
            if r_num_positions == 0:
                logger.warning(
                    f"No R-group attachment location founds for scaffold {scaffold}. "
                    + "Ensure scaffold R groups are labeled like [*:1], [*:2], etc..., "
                    + "not [*1], [*2]"
                )
                failed_attempts += 1
                scaffold_id = (scaffold_id + 1) % len(scaffolds)
                continue
            random_r_groups = self.rng.choice(
                r_groups, size=r_num_positions, replace=True
            )
            mol = utils.create_molecule(scaffold, random_r_groups)
            if mol.smiles not in generated_molecules:
                generated_molecules.add(mol.smiles)
                molecules.append(mol)
            else:
                failed_attempts += 1
            scaffold_id = (scaffold_id + 1) % len(scaffolds)
        if failed_attempts == 1000:
            logger.warning(f"Could not generate {num_mols} for some reason.")
        if len(molecules) == 0:
            raise ValueError(
                "Alert: No molecules were produced during random generation. Please"
                + " ensure the scaffolds and R-groups are correctly implemented."
            )
        if self.verbose > 0:
            print()

        return molecules

    def mate(
        self, r_groups: list[str], scaffolds: list[str], num_mols: int
    ) -> list[Molecule]:
        """Create a list of pydantic containing num_mols molecule
        chosen randomly from the r-groups and scaffolds

        Checks for any errors in R-group addition and then
        logger.debugs a statement declaring if it found one

        Args:
            r_groups (list[str]):
                list of all R-groups
            scaffolds (list[str]):
                list of all scaffolds
            num_mols (int):
                The number of Molecules to create

        Returns:
            chem_list(list[Molecule]):
                a list containing the Molecule values.
        """

        molecules = []  # stores all the newly generated molecules
        generated_molecules = set()  # ensure we only create a new molecule once
        unique_scaffolds = set(scaffolds)
        r_index_spots = 0
        for sub_i in unique_scaffolds:
            r_index_spots = r_index_spots + sub_i.count("[*:")
        total_diff_comb = len(r_groups) ** r_index_spots
        if total_diff_comb < num_mols:
            logger.info(f"changing num children from {num_mols} to {total_diff_comb}")
            num_mols = total_diff_comb

        counter = 0
        while len(generated_molecules) < num_mols:
            scaffold_index = self.rng.integers(0, len(scaffolds))
            sub = scaffolds[scaffold_index]
            r_num_positions = sub.count("[*:")
            random_r_list = [
                self.rng.integers(0, len(r_groups)) for _ in range(r_num_positions)
            ]
            random_tuple = [scaffold_index, random_r_list]
            final_chem = sub
            for r_pos_index, r_index in enumerate(random_r_list, start=1):
                R = r_groups[r_index]
                final_chem = utils.insert_r_group(final_chem, r_pos_index, R)

            if final_chem not in generated_molecules:
                final_chem_mol = Chem.MolFromSmiles(final_chem)
                if final_chem_mol is not None:
                    generated_molecules.add(final_chem)
                    r_groups_dict = {
                        f"R{i+1}": r_groups[random_tuple[1][i]]
                        for i in range(r_num_positions)
                    }
                    molecules.append(
                        Molecule(
                            smiles=final_chem,
                            structure_scaffold=sub,
                            r_groups=r_groups_dict,
                        )
                    )
                else:
                    logger.info(f"R-group could not be inserted: {sub} and R-group {R}")

            counter = counter + 1
            if counter > (num_mols * 5):
                num_mols = num_mols - 1
                counter = 0
        return molecules

    def reproduce(self, parent_list: list[Polymer], num_children_edit: int = -1):
        """Creates the next generation of molecules/polymers.

        Crosses over the parents, and mutates the children.

        Args:
            parent_list (list: Polymer):
                a list containing the dictionaries of the parent molecules

            num_children_edit (int):
                number of children per family, in case that needed to be changed.
                self.num_children_per_family is default

        Returns:
            all_children_mutated (list: Molecule):
                list of Molecules of the next generation (children) of the parents
        """
        if num_children_edit > 0:
            self.num_children_per_family = num_children_edit
        all_children = self.crossover_all(parent_list)
        all_children = self.mutate_all(all_children, fraction_mutation=0.2)
        # some indices can be None here
        children_poly = self.polymerize_all(all_children)
        indices_to_remove = set()
        max_attempts = 3
        for i in range(len(children_poly)):
            pol = children_poly[i]
            mol = all_children[i]
            attempts = 0
            while (
                pol is None or pol.canonical_smiles in self.polymer_dict
            ) and attempts < max_attempts:
                mol = self.mutate(mol)
                pol = self.polymerize_with_random_rxn_procedure(mol)
                attempts += 1
            if attempts >= max_attempts or pol is None:
                indices_to_remove.add(i)
            else:
                children_poly[i] = pol
                all_children[i] = mol

        children_poly = [
            pol for i, pol in enumerate(children_poly) if i not in indices_to_remove
        ]
        all_children = [
            mol for i, mol in enumerate(all_children) if i not in indices_to_remove
        ]

        store_poly_chem.store_molecules(self, all_children)
        self.current_population = all_children
        self.generation = self.generation + 1
        if len(self.current_population) == 0:
            raise ValueError(
                "No children exist in the current generation. This means all "
                + "generated children either failed to polymerize or had been seen "
                + "in a prior generation. Try adding more R-groups or templates."
            )
        return children_poly

    def crossover_all(self, molecule_list: list[Polymer]):
        """Iterates through the parents and mates them based on the
         partner_selection variable

        Args:
            molecule_list (list: Polymer):
                list containing the dictionaries of the parent molecules to be mated

        Return:
            total_children_list (list: Molecule):
               list containing all the children after crossing over the

        Raises:
            ValueError: If the partner_selection value is invalid.
        """
        total_children_list = []
        total_child_smiles = []

        children_num = self.num_children_per_family
        if self.partner_selection == "all":
            for i, chem in enumerate(molecule_list):
                for j in range(i + 1, len(molecule_list)):
                    children_list = self.crossover(
                        molecule_list[i].molecule,
                        molecule_list[j].molecule,
                        children_num,
                    )
                    total_children_list.extend(children_list)

        elif self.partner_selection == "random":
            total_parent_array = []
            prog = 0
            for index_family in range(0, self.num_families):
                if self.verbose > 0:
                    PROG = 100 * round(prog / self.num_families, 3)
                    sys.stdout.write("Crossing-Over, \r%d%% done" % PROG)
                    sys.stdout.flush()
                    prog += 1
                # logger.debug ("test in crossover")
                parent_array = self.rng.choice(molecule_list, size=2, replace=False)
                parent_array = parent_array.tolist()
                if parent_array not in total_parent_array:
                    children_list = self.crossover(
                        parent_array[0].molecule, parent_array[1].molecule, children_num
                    )

                    for child in children_list:
                        if child.smiles not in total_child_smiles:
                            total_child_smiles.append(child.smiles)
                            total_children_list.append(child)

                    total_parent_array.append(parent_array)

        else:
            raise ValueError(
                "Please choose a valid selection scheme. {} invalid.".format(
                    self.partner_selection
                )
            )
        return total_children_list

    def crossover(self, parent_1: Molecule, parent_2: Molecule, children_num: int):
        """Takes two parents and mates them.
        considers scaffolds as gender, and then randomly chooses the R-groups

        Args:
            parent_1 (Molecule):
                Molecules of parent 1 to be mated with parent_2

            parent_2 (Molecule):
                Molecules of parent 2 to be mated with parent_1

            children_num (int):
                number of children to be reproduced for this family

        Returns:
            list_children (list: Molecule):
                list of the children of this family
        """
        list_children = []
        tot_r_group = []
        r_num_positions_1 = parent_1.structure_scaffold.count("[*:")
        r_num_positions_2 = parent_2.structure_scaffold.count("[*:")
        for r_pos_index in range(1, r_num_positions_1 + 1):
            r_string = "R" + str(r_pos_index)
            tot_r_group.append(parent_1.r_groups[r_string])
        for r_pos_index in range(1, r_num_positions_2 + 1):
            r_string = "R" + str(r_pos_index)
            tot_r_group.append(parent_2.r_groups[r_string])

        tot_r_group = list(set(tot_r_group))
        list_children = self.mate(
            tot_r_group,
            [parent_1.structure_scaffold, parent_2.structure_scaffold],
            children_num,
        )
        return list_children

    def mutate_all(
        self, molecule_list: list[Molecule], fraction_mutation: float = 0.075
    ):
        """mutates a certain percenatge (normalky distrbiuted) of the entire children
        list (at least 1 and most all)

        Args:
            molecule_list (list: Molecule):
                The total children list where some will be mutated

            fraction_mutations (float):
                The mean value of the normal distribution

        Returns:
            Molecule_list  (list: Molecule):
                New list of children where some have been mutated
        """
        mutation_sigma = 20
        num_to_mutate = min(
            max(
                round(
                    self.rng.normal(
                        len(molecule_list) * fraction_mutation, mutation_sigma
                    )
                ),
                1,
            ),
            len(molecule_list),
        )

        num_to_mutate = abs(num_to_mutate)
        num_to_mutate = min(len(molecule_list), num_to_mutate)
        indices = [x for x in range(len(molecule_list))]
        molecules_to_mutate = list(
            self.rng.choice(indices, size=num_to_mutate, replace=False)
        )
        for mol_ndex in molecules_to_mutate:
            molecule_list[mol_ndex] = self.mutate(molecule_list[mol_ndex])
        return molecule_list

    def mutate(self, molecule: Molecule):
        """Mutations can be:
        1. change in R-group (drop a letter, add a new R group from the list)
        2. Change in the position of R-group

        Mutates a single molecule
        Args:
            molecule (dict):
                the molecule to be mutated

        Returns:
            molecule (dict):
                the molecule after being mutated
        """
        old_smiles = molecule.smiles
        for _ in range(30):
            r_num_positions = molecule.structure_scaffold.count("[*:")
            R_pos_to_mutate = self.rng.choice(range(1, r_num_positions + 1))
            R_to_mutate = self.rng.choice(self.r_groups)
            molecule = self.exchange(
                molecule, r_pos=R_pos_to_mutate, new_r_group=R_to_mutate
            )
            if old_smiles != molecule.smiles:
                break
        return molecule

    def exchange(self, molecule: Molecule, r_pos: int = 1, new_r_group: str = "[*][H]"):
        """Swaps the R-group with another inputted one at a specific location

        Args:
            molecule (Molecule):
                the dictionary of the molecule to be altered

            r_pos (int):
                the position of R_group to be replaced

            new_r_group (str):
                The new R-group to replace the old one

        Returns:
            molecule (Molecule):
                return this if swap completed successfully

            molecule_old (Molecule):
                otherwise returns the old molecule
        """
        molecule_old = molecule.copy()
        temp_r_group_dict = molecule.r_groups.copy()
        r_string_index = "R" + str(r_pos)
        temp_r_group_dict[r_string_index] = new_r_group
        scaffold = molecule.structure_scaffold
        r_num_positions = scaffold.count("[*:")
        mol = scaffold
        for r_pos_index in range(1, r_num_positions + 1):
            r_string_index = "R" + str(r_pos_index)
            R = temp_r_group_dict[r_string_index]
            mol = utils.insert_r_group(mol, r_pos_index, R)
            if mol is None:
                break
        if mol is not None:
            molecule.smiles = mol
            r_string_index = "R" + str(r_pos)
            molecule.r_groups[r_string_index] = new_r_group
            return molecule
        else:
            logger.info(
                "In mutation --> R-group could not be inserted: "
                + f"{scaffold},R-group: {R}, returning old molecules "
            )
            return molecule_old

    def polymerize_all(self, molecule_list: list[Molecule]) -> list[Polymer]:
        """Polymerizes and fingerprints all the molecules with the R-groups in the
        dictionary, and saves the "polymers" in image directory.

        Checks for any errors in polymerization and then
        logger.debugs a statement declaring if it found one

        Args:
            molecule_list (list:Molecule):
                a list containing the Molecule values.

        Returns:
            (list:Polymer) :
                list containing the Polymers of the molecules along with the
                fingerprints.
        """
        polymer_list = []
        prog = 0
        for mol in molecule_list:
            if self.verbose > 0:
                PROG = 100 * round(prog / (len(molecule_list)), 3)
                sys.stdout.write("Polymerization: \r%d%% done" % PROG)
                sys.stdout.flush()
                prog += 1
            pol = self.polymerize_with_random_rxn_procedure(mol)
            # still add polymer to list, will have to sort out later though
            polymer_list.append(pol)
            if pol is None:
                logger.info(
                    "Substructure could not be polymerized with these R-group: "
                    + f"{mol.structure_scaffold} -> {mol.smiles}"
                )
        return polymer_list

    def fingerprint_efficient(self, polymer_list) -> list[Polymer]:
        """Makes sure that only new polymers are fingerprinted

        Checks for any errors in polymerization and then
        logger.debugs a statement declaring if it found one

        Args:
            polymer_list (list:polymer):
                a list containing the polymer values.

        Returns:
            (list:Polymer) :
                list containing the Polymers of the molecules along with the
                fingerprints.
        """
        polymer_to_fingerprint = []
        already_fingerprint = []
        for poly in polymer_list:
            if poly.canonical_smiles not in self.polymer_fingerprint:
                polymer_to_fingerprint.append(poly)
            else:
                poly.fingerprint = self.polymer_fingerprint[poly.canonical_smiles]
                already_fingerprint.append(poly)
        polymer_to_fingerprint = self.fp_function(polymer_list)
        polymer_list = already_fingerprint + polymer_to_fingerprint
        return polymer_list

    def polymerize_with_random_rxn_procedure(self, molecule: Molecule) -> Polymer:
        """Polymerize molecule based on random reaction procedure from scaffold dict"""
        rxn_procedure_name = self.rng.choice(
            self.scaffold_reaction_dict[molecule.structure_scaffold]
        )
        rxn_procedure = self.reaction_procedure_dict[rxn_procedure_name]
        return self.polymerize(molecule, rxn_procedure)

    def polymerize(
        self, molecule: Molecule, rxn_procedure: ReactionProcedure
    ) -> Polymer:
        """Runs a polymerization reaction on a molecule.

        Runs the reaction scaffold, and then uses that to polymerize an
        inputted molecule

        Args:
            molecule (Molecule):
                the molecule to be polymerized
        Returns:
            Polymer if success, else None
        """
        products = Chem.MolFromSmiles(molecule.smiles)
        for step in rxn_procedure.steps:
            rxn = rdChemReactions.ReactionFromSmarts(step)
            reactants = [products]
            # get additional reactants needed in the template. First one is the
            # original molecule
            for ind in range((rxn.GetNumReactantTemplates()) - 1):
                reactants.append(rxn.GetReactantTemplate(ind + 1))
            try:
                products = rxn.RunReactants(tuple(reactants))
                products = products[0][0]
                Chem.SanitizeMol(products)
                products.UpdatePropertyCache()
            except Exception as e:
                logger.error(f"Error during polymerization of {molecule.smiles}: {e}")
                return None

        if products is None:
            logger.warning(f"Failed to polymerize {molecule.smiles}")
            return None
        smiles = Chem.MolToSmiles(products)
        try:
            canonical_smiles = canonicalize(smiles)
        # certain RDKit exceptions can be thrown that fail this
        except Exception as e:
            logger.warning(f"Failed canonicalizing {smiles}: {e}")
            return None
        # if canonical_smiles not in or in self.polymer_fingerprint:
        return Polymer(
            molecule=molecule,
            smiles=smiles,
            canonical_smiles=canonical_smiles,
            reaction_proc=rxn_procedure,
        )
