import logging, sys
import pymol2
from rdkit import Chem
import math
from tqdm import tqdm
import matplotlib
import matplotlib.font_manager

matplotlib.font_manager._load_fontmanager(try_read_cache=True)
matplotlib.set_loglevel("critical")

from choppa.utils import (
    show_contacts,
    get_ligand_resnames_from_pdb_str,
    split_pdb_str,
    get_contacts_mda,
    biopython_to_mda,
)
from choppa.logoplots import (
    LogoPlot,
    WHITE_EMPTY_SQUARE,
    render_singleres_logoplot,
)
from choppa.data.templates.resources import HTML_TEMPLATE, LOGOPLOT_TEMPLATE
from jinja2 import Template


logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()

HEX_COLOR_CODES = [
    "#ffffff",
    "#ff9e83",
    "#ff8a6c",
    "#ff7454",
    "#ff5c3d",
    "#ff3f25",
    "#ff0707",
    "#642df0",
]


class PublicationView:
    """
    Uses the PyMOL API to create a session file for publication-ready views of the fitness data on
    top of the complex PDB. Users will need to `ray` the PyMOL view with the desired `ray`
    settings in the GUI application themselves, but a combination of pre-set `ray` settings is provided
    in the `.pse` file.
    """

    def __init__(
        self,
        filled_aligned_fitness_dict,
        complex,
        complex_rdkit,
        fitness_threshold,
        output_session_file="out.pse",
    ):
        self.fitness_dict = filled_aligned_fitness_dict
        self.complex = complex
        self.fitness_threshold = fitness_threshold
        self.output_session_file = output_session_file
        self.ligand_present = False  # will turn to True if found.

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
        p.cmd.read_pdbstr(self.complex_pdb_str, "complex")

        if remove_solvents:
            p.cmd.remove("inorganic")
            p.cmd.remove("solvent")

    def count_fit_residues(self, fitness_data):
        """
        For a dict with mutants in a fitness dict, counts the number of fit mutants
        """
        fit_mutants_count = 0
        wildtype = fitness_data["wildtype"]["aa"]
        for mut_dict in fitness_data["mutants"]:
            if mut_dict["aa"] == wildtype:
                continue  # skip wildtype, don't want to count this as a mutant

            if mut_dict["fitness"] >= self.fitness_threshold:
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
            if "mutants" in fitness_data:
                # set the number of counted fit mutants to a string variable that we'll
                # parse in the other dict
                residue_mutability_levels[
                    f"n_fit_{self.count_fit_residues(fitness_data)}"
                ].append(str(i))
            else:
                residue_mutability_levels["no_fitness_data"].append(str(i))

        # make the index values separated by '+' so that PyMOL can read it
        for k, v in residue_mutability_levels.items():
            residue_mutability_levels[k] = "+".join(v)
        logger.info(
            f"PyMOL session: fitness degree per residue found using threshold {self.fitness_threshold}:\n{residue_mutability_levels}\n"
        )

        mutability_color_dict = (
            {  # TODO: convert color coding to match 3DMol color specs
                "n_fit_0": "white",
                "n_fit_1": "salmon",
                "n_fit_2": "deepsalmon",
                "n_fit_3": "tv_red",
                "n_fit_4": "red",
                "n_fit_5": "firebrick",
                "no_fitness_data": "skyblue",
            }
        )
        return residue_mutability_levels, mutability_color_dict

    def pymol_color_by_fitness(self, p):
        """
        With a pymol session set up with a system using `self.pymol_setup_system()`,
        integrates `fitness` data by coloring residues by mutability degree.
        """
        logger.info("PyMOL session: coloring system surface with fitness data..")

        # first get the fitness degree per residue and what colors to make them
        residue_mutability_levels, mutability_color_dict = self.pymol_color_coder()
        self.residue_mutability_levels = residue_mutability_levels

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
            self.ligand_present = True
            # need to use the PyMOL DSL for selection. Construct the string first.
            ligand_selector = f"resn {ligands[0]}"
            for lig in ligands[1:]:  # add additional ligand entries, if there are any
                ligand_selector += f" and resn {lig}"
        else:
            logger.warn("WARNING: no ligands detected in PDB file.")

        # make the selections
        if self.ligand_present:
            p.cmd.select("ligand", ligand_selector)
        p.cmd.select(  # PyMOL will have made only the protein cartoon, so can just select that way
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
        p.cmd.set("ray_trace_mode", 1)  # smooth with black outlining on objects
        p.cmd.set(
            "specular", 0
        )  # removes reflections on surface for easier interpretation

        if ligands_in_system:
            # select the ligand and subpocket residues, show them as sticks w/o nonpolar Hs
            # TODO: set ligand stick color that conforms to HTML view
            p.cmd.util.cba(144, "ligand", _self=p.cmd)
            p.cmd.show("sticks", "ligand")
            p.cmd.show("spheres", "ligand")
            p.cmd.set("stick_radius", "0.15")
            p.cmd.set("sphere_scale", "0.15")

    def pymol_add_interactions(self, p, mutability_color_dict):
        """
        Adds interactions to pymol session if a ligand is present. Interactions are colored by
        fitness of contacted residues, not by interaction type.
        """
        logger.info(
            "PyMOL session: adding ligand-protein interactions (contacts) colored by fitness degree"
        )

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


class InteractiveView:
    """
    Uses 3DMol and Jinja to create a single HTML file that can be hosted anywhere to enable shareable
    interactive views of the fitness data on top of the complex PDB
    """

    def __init__(
        self,
        filled_aligned_fitness_dict,
        complex,
        complex_rdkit,
        fitness_threshold,
        output_session_file="out.html",
    ):
        self.fitness_dict = filled_aligned_fitness_dict
        self.complex = complex
        self.fitness_threshold = fitness_threshold
        self.output_session_file = output_session_file
        self.ligand_present = False  # will turn to True if found.

        # get the PDB file as a string from RDKit
        self.complex_pdb_str = Chem.MolToPDBBlock(complex_rdkit)

    def get_confidence_limits(self):
        """
        Figures out what the global maximum and minimum is of the confidence measures (e.g. number of reads)
        in the experimental protocol of the fitness data. If there is no confidence measure, returns False
        """
        self.confidence = False
        confidence_values = []

        for i, res in self.fitness_dict.items():
            if "mutants" in res:  # skips over PDB residues that don't have fitness data
                mut_conf_values = [mut["confidence"] for mut in res["mutants"]]
                wildtype_conf_value = res["wildtype"]["confidence"]

                for conf_val in mut_conf_values + [wildtype_conf_value]:
                    confidence_values.append(conf_val)
        if math.isnan(confidence_values[0]):
            return False, False
        else:
            return [min(confidence_values), max(confidence_values)]

    def get_logoplot_dict(self, confidence_lims, multiprocess=True, max_workers=None):
        """
        For a fitness dict, load all base64 logoplots into memory using multithreading if requested.

        Instead of adding base64 strings to fitness dict (making it uninterpretable), make a separate
        dict that mimics the form of fitness dict.
        """
        logger.info(f"Generating logoplots for {len(self.fitness_dict)} residues.")
        logoplot_dict = {}

        if multiprocess:
            logger.info("Using MultiProcessing")
            from concurrent.futures import ProcessPoolExecutor

            # use multithreading to generate logoplots
            with ProcessPoolExecutor(
                max_workers=max_workers
            ) as executor:  # Use optimal thread pool size
                futures = {
                    executor.submit(
                        self._make_logoplot_residue,
                        idx,
                        residue_fitness_dict,
                        confidence_lims,
                    )
                    for idx, residue_fitness_dict in self.fitness_dict.items()
                }

            with tqdm(total=len(futures)) as progress_bar:
                for future in futures:
                    try:
                        logoplot, idx = future.result()
                        logoplot_dict[idx] = logoplot
                        progress_bar.update(
                            1
                        )  # Update progress bar for each completed task
                    except Exception as e:
                        # Handle potential exceptions raised during thread execution
                        logger.error(f"Error generating logoplot: {e}")
                        progress_bar.update(1)  # Update progress bar even on error

        else:
            for idx, residue_fitness_dict in tqdm(self.fitness_dict.items()):
                logoplot, idx = self._make_logoplot_residue(
                    idx, residue_fitness_dict, confidence_lims
                )
                logoplot_dict[idx] = logoplot

        # sort the logoplot dict by residue index so that the multiprocessed logoplots are in same order as
        # if they were generated sequentially.
        logoplot_dict = dict(sorted(logoplot_dict.items()))

        return logoplot_dict

    def _make_logoplot_residue(self, idx, residue_fitness_dict, confidence_lims):
        """
        For a single residue, generate logoplots for wildtype, fit and unfit mutants.
        """
        if not "aa" in residue_fitness_dict["wildtype"]:
            return (
                {
                    "fitness_aligned_index": idx,
                    "fitness_csv_index": idx,
                    "logoplots_base64": {
                        "wildtype": render_singleres_logoplot(
                            "".join(residue_fitness_dict["wildtype"])
                        ),
                        "fit": WHITE_EMPTY_SQUARE,
                        "unfit": WHITE_EMPTY_SQUARE,
                    },
                },
                idx,
            )

        wildtype_base64, fit_base64, unfit_base64 = LogoPlot(
            residue_fitness_dict, fitness_threshold=self.fitness_threshold
        ).build_logoplot(
            global_min_confidence=confidence_lims[0],
            global_max_confidence=confidence_lims[1],
        )

        return (
            {
                "fitness_aligned_index": residue_fitness_dict["fitness_aligned_index"],
                "fitness_csv_index": residue_fitness_dict["fitness_csv_index"],
                "logoplots_base64": {
                    "wildtype": wildtype_base64,
                    "fit": fit_base64,
                    "unfit": unfit_base64,
                },
            },
            idx,
        )

    def get_surface_coloring_dict(self):
        """
        Based on fitness coloring, creates a dict where keys are colors, values are residue numbers.
        """
        color_res_dict = {color: [] for color in HEX_COLOR_CODES}
        for resindex, fitness_data in self.fitness_dict.items():
            if not "mutants" in fitness_data:
                # no fitness data for this residue after alignment between Fitness CSV and input PDB
                color_res_dict["#642df0"].append(resindex)  # makes residue surface blue
                continue
            fit_mutants = [
                mut
                for mut in fitness_data["mutants"]
                if mut["fitness"] > self.fitness_threshold
            ]
            if (
                len(fit_mutants) <= 6
            ):  # color by increasing redness the more fit mutants there are
                color_res_dict[HEX_COLOR_CODES[len(fit_mutants)]].append(resindex)
            else:  # if there are more than 5 fit mutants just color it the most red - super mutable in any case
                color_res_dict[HEX_COLOR_CODES[6]].append(resindex)

        return color_res_dict

    def surface_coloring_dict_to_js(self, color_res_dict):
        """
        Transforms a dictionary of residue indices per color (hex) to a JavaScript-compatible
        string.
        """
        residue_coloring_function_js = ""
        start = True
        for color, residues in color_res_dict.items():
            residues = [
                f"'{res}'" for res in residues
            ]  # need to wrap the string in quotes *within* the JS code
            if start:
                residue_coloring_function_js += (
                    "if (["
                    + ",".join(residues)
                    + "].includes(atom_residx)){ \n return '"
                    + color
                    + "' \n "
                )
                start = False
            else:
                residue_coloring_function_js += (
                    "} else if (["
                    + ",".join(residues)
                    + "].includes(atom_residx)){ \n return '"
                    + color
                    + "' \n "
                )
        return residue_coloring_function_js

    def get_interaction_dict(self):
        """
        Generates interactions to be displayed on the interactive HTML view. Interactions are colored
        by the same rules as for PyMOL (`render.PublicationView()`), but a dict of interactions is
        used which is generated in `render.PublicationView().pymol_add_interactions()`.
        """
        intn_dict = {}
        intn_count = 0
        get_contacts_mda(self.complex)

        for contact in get_contacts_mda(self.complex):
            # first find the color.
            contact_color = [
                k
                for k, v in self.get_surface_coloring_dict().items()
                if contact[1].resid in v
            ]

            # check if the residue atom is in backbone, if so just overwrite the color to make the
            # contact color green.
            backbone_atoms = [
                res
                for res in biopython_to_mda(self.complex).select_atoms(
                    "protein and backbone"
                )
            ]
            if contact[1] in backbone_atoms:
                contact_color = ["#047b14"]  # green

            # then get the coordinates for the atoms involved; add to intn_dict
            # while also prepending with an index so that we don't overwrite
            # interactions purely keyed on residue alone (one res can have multiple
            # interactions).
            intn_dict[f"{intn_count}_{contact[1].resname}{contact[1].resid}"] = {
                "lig_at_x": contact[0].position[0],
                "lig_at_y": contact[0].position[1],
                "lig_at_z": contact[0].position[2],
                "prot_at_x": contact[1].position[0],
                "prot_at_y": contact[1].position[1],
                "prot_at_z": contact[1].position[2],
                "color": contact_color[0],
            }
            intn_count += 1

        return intn_dict

    def inject_stuff_in_template(
        self, sdf_str, pdb_str, surface_coloring, logoplot_dict
    ):
        """ "
        Replaces parts of a template HTML with relevant bits of data to get to a HTML view
        of the (ligand-) protein, its fitness and its interactions (if any).

        Uses Jinja2 templating to render based on static HTML template.
        """
        # create a bunch of DIVs of the logoplots.
        logoplot_divs = ""
        for _, logoplot_data in logoplot_dict.items():
            wt = logoplot_data["logoplots_base64"]["wildtype"]
            fit = logoplot_data["logoplots_base64"]["fit"]
            unfit = logoplot_data["logoplots_base64"]["unfit"]
            data = {
                "ID": logoplot_data["fitness_aligned_index"],
                "WT_BASE64_INSERT": wt.decode("ascii"),
                "FIT_BASE64_INSERT": fit.decode("ascii"),
                "UNFIT_BASE64_INSERT": unfit.decode("ascii"),
            }
            logoplot_template_str = open(LOGOPLOT_TEMPLATE).read()
            logoplot_template = Template(logoplot_template_str)
            logoplot_divs += logoplot_template.render(data)

        # add the PDB (protein) and SDF (ligand)
        intr_dct = str(self.get_interaction_dict()) if self.ligand_present else "{}"

        data = {
            "PDB_INSERT": pdb_str,
            "SDF_INSERT": sdf_str,
            "LOGOPLOTS_INSERTS": logoplot_divs,
            "SURFACE_COLOR_INSERT": surface_coloring,
            "INTN_DICT_INSERT": intr_dct,
        }
        # render the template with Jinja
        template_str = open(HTML_TEMPLATE).read()
        template = Template(template_str)
        with open(self.output_session_file, "wt") as fout:
            fout.write(template.render(data))

    def render(self):
        # check if we have confidences, if we do then record the [min, max]
        confidence_lims = self.get_confidence_limits()

        # load logoplots into memory as base64 images in a dict. This will take 3-10MB of memory,
        # but may be more if the PDB is huge. We could opt to do this later but we'll have to
        # load it into memory at some point anyway to be able to write out the final HTML.
        # Plus, this makes the multithreading easier to debug as there are fewer layers.
        logoplot_dict = self.get_logoplot_dict(confidence_lims)

        # # get the strings for the PDB (prot) and the SDF (lig, if present)
        lig_sdf_str, prot_pdb_str = split_pdb_str(self.complex_pdb_str)
        if lig_sdf_str:
            self.ligand_present = True
        # # define the coloring of the surface
        surface_coloring = self.surface_coloring_dict_to_js(
            self.get_surface_coloring_dict()
        )

        # do a dirty HTML generation using the logoplot and fitness dicts.
        self.inject_stuff_in_template(
            lig_sdf_str, prot_pdb_str, surface_coloring, logoplot_dict
        )


if __name__ == "__main__":
    from data.toy_data.resources import (
        TOY_COMPLEX,
        TOY_FITNESS_DATA_COMPLETE,
        TOY_FITNESS_DATA_TRUNCATED,
    )
    from data.toy_data.resources import (
        TOY_FITNESS_DATA_SECTIONED,
        TOY_FITNESS_DATA_COMPLETE_NOCONF,
    )

    from IO.input import FitnessFactory, ComplexFactory

    fitness_dict = FitnessFactory(
        TOY_FITNESS_DATA_TRUNCATED, confidence_colname="confidence"
    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    from align import AlignFactory

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    PublicationView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()

    InteractiveView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()
