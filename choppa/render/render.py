import logging, sys
import pymol2
from choppa.render.utils import show_contacts

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()

class PYMOL():
    """
    Uses the PyMOL API to create a session file for publication-ready views of the fitness data on 
    top of the complex PDB. Users will need to `ray` the PyMOL view with the desired `ray`
    settings in the GUI application themselves, but a combination of pre-set `ray` settings is provided
    in the `.pse` file. 
    """
    def __init__(self, filled_aligned_fitness_dict, complex, complex_pdb_str, fitness_threshold):
        self.fitness_dict = filled_aligned_fitness_dict
        self.complex = complex
        self.complex_pdb_str = complex_pdb_str
        self.fitness_threshold = fitness_threshold

    def pymol_start_session(self):
        """
        Boots up a session with the publicly available PyMOL API.
        """
        p = pymol2.PyMOL()
        p.start()
        return p
    
    def pymol_setup_system(self, p, remove_solvents=True):
        """
        Sets up the protein (and ligands if present) without changing any of the looks. Also removes
        some stuff we're not interested in.
        """
        p.cmd.read_pdbstr(self.complex_pdb_str, 'complex')

        if remove_solvents:
            p.cmd.remove("inorganic")
            p.cmd.remove("solvent")

    def count_fit_residues(self, fitness_data):
        """
        For a dict with mutants in a fitness dict, counts the number of fit mutants
        """
        fit_mutants_count = 0
        wildtype = fitness_data['wildtype']['aa']
        for mut_dict in  fitness_data['mutants']:
            if mut_dict['aa'] == wildtype:
                continue # skip wildtype, don't want to count this as a mutant
            
            if mut_dict['fitness'] >= self.fitness_threshold:
                fit_mutants_count += 1

        # return while clipping to 5. At 5 fit mutants the residue is so mutable that we don't need to
        # make it any more red, it's already super bad news.
        return min(5, fit_mutants_count)


    def pymol_color_coder(self):
        """
        Given the aligned fitness dict, returns two dicts:
        - residue_mutability_levels : `{residue index : number of fit mutations, ..}`
        - mutability_color_dict : `{number of fit mutations : color, ..}`
        """

        # instantiate a dict to populate.
        residue_mutability_levels = {}
        for k in [
            "n_fit_0",
            "n_fit_1",
            "n_fit_2",
            "n_fit_3",
            "n_fit_4",
            "n_fit_5",
            "no_fitness_data",
        ]:
            residue_mutability_levels[k] = [] 

        for i, fitness_data in self.fitness_dict.items():
            if 'mutants' in fitness_data:
                # set the number of counted fit mutants to a string variable that we'll 
                # parse in the other dict
                residue_mutability_levels[f"n_fit_{self.count_fit_residues(fitness_data)}"].append(str(i))
            else:
                residue_mutability_levels["no_fitness_data"].append(str(i))
        
        # make the index values separated by '+' so that PyMOL can read it
        for k,v in residue_mutability_levels.items():
            residue_mutability_levels[k] = "+".join(v)

        mutability_color_dict = { # TODO: convert color coding to match 3DMol color specs
            "n_fit_0" : "white",
            "n_fit_1" : "salmon",
            "n_fit_2" : "deepsalmon",
            "n_fit_3" : "tv_red",
            "n_fit_4" : "red",
            "n_fit_5" : "firebrick",
            "no_fitness_data" : "gray",
        }
        return residue_mutability_levels, mutability_color_dict

    def pymol_color_by_fitness(self, p):
        """
        With a pymol session set up with a system using `self.pymol_setup_system()`,
        integrates `fitness` data by coloring residues by mutability degree.
        """
        # first get the fitness degree per residue and what colors to make them
        residue_mutability_levels, mutability_color_dict = self.pymol_color_coder()

        # now select the residues, name them and color them.
        for fitness_degree, residues in residue_mutability_levels.items():
            if residues:
                # found fitness data (i.e. residue indices) for this degree of mutability
                p.cmd.select(
                    fitness_degree,
                    f"complex and resi {residues} and polymer.protein",
                )
            else:
                # no fitness data for this level of mutability - remove from dict so 
                # that PyMOL doesn't select the residue.
                mutability_color_dict.pop(fitness_degree)

        for fitness_degree_name, color in mutability_color_dict.items():
            p.cmd.set("surface_color", color, f"({fitness_degree_name})")

    def pymol_prettify_system(self, p):
        """
        With a pymol session set up with a system using `self.pymol_setup_system()`,
        makes the session pretty.
        """
                
        # tmp
        p.cmd.show("surface")

    def pymol_add_interactions(self, p):
        """
        Adds interactions to pymol session if a ligand is present. Interactions are colored by
        fitness of contacted residues, not by interaction type.
        """
    
    def pymol_write_session(self, p, out_filename):
        """
        Writes out a pymol session to a `.pse` file.
        """
        p.cmd.save(out_filename)

    def render(self):
        """
        Renders the PyMOL session file.
        """
        # start the session
        p = self.pymol_start_session()

        # set up the protein(-ligand) into the session
        self.pymol_setup_system(p)

        # color it by fitness
        self.pymol_color_by_fitness(p)

        # make the view pretty
        self.pymol_prettify_system(p)

        # finally write to a session file (`.pse`)
        self.pymol_write_session(p, "test_sess.pse")

class HTML():
    """
    Uses 3DMol and Jinja to create a single HTML file that can be hosted anywhere to enable shareable 
    interactive views of the fitness data on top of the complex PDB
    """


if __name__ == "__main__":
    from choppa.data.toy_data.resources import TOY_COMPLEX, TOY_FITNESS_DATA_COMPLETE, TOY_FITNESS_DATA_TRUNCATED
    from choppa.data.toy_data.resources import TOY_FITNESS_DATA_SECTIONED, TOY_FITNESS_DATA_COMPLETE_NOCONF

    from choppa.IO.input import FitnessFactory, ComplexFactory

    fitness_dict = FitnessFactory(TOY_FITNESS_DATA_SECTIONED, 
                                    # confidence_colname="confidence"
                                    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()

    from choppa.align.align import AlignFactory
    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    PYMOL(filled_aligned_fitness_dict, 
          complex, ComplexFactory(TOY_COMPLEX).load_pdb_string(), 
          fitness_threshold=0.7).render()