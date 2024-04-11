import logging, sys
import pymol2
from io import StringIO
from rdkit import Chem

from choppa.render.utils import show_contacts, get_ligand_resnames_from_pdb_str

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()

class PYMOL():
    """
    Uses the PyMOL API to create a session file for publication-ready views of the fitness data on 
    top of the complex PDB. Users will need to `ray` the PyMOL view with the desired `ray`
    settings in the GUI application themselves, but a combination of pre-set `ray` settings is provided
    in the `.pse` file. 
    """
    def __init__(self, filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold, output_session_file="out.pse"):
        self.fitness_dict = filled_aligned_fitness_dict
        self.complex = complex
        self.fitness_threshold = fitness_threshold
        self.output_session_file = output_session_file

        # get the PDB file as a string from RDKit
        self.complex_pdb_str = Chem.MolToPDBBlock(complex_rdkit)

    def pymol_start_session(self):
        """
        Boots up a session with the publicly available PyMOL API.
        """
        logger.info("Starting PyMOL session")
        p = pymol2.PyMOL()
        p.start()
        return p
    
    def pymol_setup_system(self, p, remove_solvents=True):
        """
        Sets up the protein (and ligands if present) without changing any of the looks. Also removes
        some stuff we're not interested in.
        """
        logger.info("PyMOL session: setting up system")
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
        logger.info(f"PyMOL session: fitness degree per residue found using threshold {self.fitness_threshold}:\n{residue_mutability_levels}\n")

        mutability_color_dict = { # TODO: convert color coding to match 3DMol color specs
            "n_fit_0" : "white",
            "n_fit_1" : "salmon",
            "n_fit_2" : "deepsalmon",
            "n_fit_3" : "tv_red",
            "n_fit_4" : "red",
            "n_fit_5" : "firebrick",
            "no_fitness_data" : "skyblue",
        }
        return residue_mutability_levels, mutability_color_dict

    def pymol_color_by_fitness(self, p):
        """
        With a pymol session set up with a system using `self.pymol_setup_system()`,
        integrates `fitness` data by coloring residues by mutability degree.
        """
        logger.info("PyMOL session: coloring system surface with fitness data..")

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
        
        return mutability_color_dict

    def pymol_select_components(self, p):
        """
        Makes selections in PyMOL for ligand, protein, binding site. Returns whether there
        is/are (a) ligand(s) present in the system.
        """
        ligands = get_ligand_resnames_from_pdb_str(self.complex_pdb_str)
        if ligands:
            # need to use the PyMOL DSL for selection. Construct the string first.
            ligand_selector = f"resn {ligands[0]}"
            for lig in ligands[1:]: # add additional ligand entries, if there are any
                ligand_selector += f" and resn {lig}"

        # make the selections
        p.cmd.select("ligand", ligand_selector)
        p.cmd.select( # PyMOL will have made only the protein cartoon, so can just select that way
            "receptor", "rep cartoon"
        )  

        return ligands
    
    def pymol_prettify_system(self, p, ligands_in_system):
        """
        With a pymol session set up with a system using `self.pymol_setup_system()`,
        makes the session pretty. This code isn't pretty though, that's just because
        of how the PyMOL API is constructed.
        """
        logger.info("PyMOL session: prettifying view")

        # reset the view
        p.cmd.select("None") 
        p.cmd.set("bg_rgb", "white")
        p.cmd.bg_color("white")
        p.cmd.hide("everything")

        # show the protein surface
        p.cmd.show("surface", "receptor")
        p.cmd.set("surface_mode", 3)
        
        # set some variables to improve the image when the user `ray`s the session once loaded
        # these increase the loading/ray time in PyMOL significantly but are needed to make 
        # publication-ready figures. These could be adjusted to conform to journals' demands.
        p.cmd.set("surface_quality", 2)
        p.cmd.set("antialias", 2)
        p.cmd.set("ray_shadows", 0) 
        p.cmd.set("ray_trace_mode", 1) # smooth with black outlining on objects
        p.cmd.set("specular", 0) # removes reflections on surface for easier interpretation

        if ligands_in_system:
            # select the ligand and subpocket residues, show them as sticks w/o nonpolar Hs
            # TODO: set ligand stick color that conforms to HTML view
            p.cmd.util.cba(144,"ligand",_self=p.cmd)
            p.cmd.show("sticks", "ligand")
            p.cmd.show("spheres", "ligand")
            p.cmd.set("stick_radius", "0.15")
            p.cmd.set("sphere_scale", "0.15")
            
    def pymol_add_interactions(self, p, mutability_color_dict):
        """
        Adds interactions to pymol session if a ligand is present. Interactions are colored by
        fitness of contacted residues, not by interaction type.
        """
        logger.info("PyMOL session: adding ligand-protein interactions (contacts) colored by fitness degree")

        # for each degree of fitness (0, 1, .. 5, no_data), show ligand-protein contacts and color them
        # in the same way that we colored the residue surface
        for fitness_degree_selection, color in mutability_color_dict.items():
            show_contacts(p, fitness_degree_selection, "ligand", contact_color=color)
            
    def pymol_write_session(self, p, out_filename):
        """
        Writes out a pymol session to a `.pse` file.
        """
        logger.info(f"PyMOL session: writing session file to {out_filename}\n")
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
        mutability_color_dict = self.pymol_color_by_fitness(p)

        # make selections and figure out whether there is a ligand present
        ligands_in_system = self.pymol_select_components(p)

        # make the view pretty
        self.pymol_prettify_system(p, ligands_in_system)
        
        if ligands_in_system:
            # add interactions
            self.pymol_add_interactions(p, mutability_color_dict)

        # finally write to a session file (`.pse`)
        self.pymol_write_session(p, self.output_session_file)

class HTML():
    """
    Uses 3DMol and Jinja to create a single HTML file that can be hosted anywhere to enable shareable 
    interactive views of the fitness data on top of the complex PDB
    """



if __name__ == "__main__":
    from choppa.data.toy_data.resources import TOY_COMPLEX, TOY_FITNESS_DATA_COMPLETE, TOY_FITNESS_DATA_TRUNCATED
    from choppa.data.toy_data.resources import TOY_FITNESS_DATA_SECTIONED, TOY_FITNESS_DATA_COMPLETE_NOCONF

    from choppa.IO.input import FitnessFactory, ComplexFactory

    fitness_dict = FitnessFactory(TOY_FITNESS_DATA_COMPLETE, 
                                    confidence_colname="confidence"
                                    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    from choppa.align.align import AlignFactory
    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    PYMOL(filled_aligned_fitness_dict, 
          complex,
          complex_rdkit,
          fitness_threshold=0.7).render()

