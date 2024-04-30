import click
import shutil
from typing import Optional

from choppa.IO.input import FitnessFactory, ComplexFactory
from choppa.align.align import AlignFactory
from choppa.render.render import PublicationView, InteractiveView
from choppa.cli.utils import SpecialHelpOrder


@click.group(
    cls=SpecialHelpOrder,
    context_settings={"max_content_width": shutil.get_terminal_size().columns - 20},
)
@click.version_option()
def cli():
    "Integrated mutational and structural biology data into a concerted HTML view."


@cli.command(
    name="render",
    help="Create fitness view as a publication-ready PyMOL session file and a read-for-sharing interactive HTML file. ",
    short_help="Create fitness view as a publication-ready PyMOL session file and a read-for-sharing interactive HTML file.",
)
@click.option(
    "-p",
    "--pdb-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, writable=True),
    help="Path to a PDB file to create fitness view for.",
    required=True,
)
@click.option(
    "-f",
    "--fitness-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, writable=True),
    help="Path to a CSV file with fitness data to create fitness view for.",
    required=True,
)
@click.option(
    "-ft",
    "--fitness-threshold",
    type=click.FLOAT,
    help="Fitness threshold to determine whether a mutant is fit or not.",
    required=True,
)
@click.option(
    "-op",
    "--outfile-publication",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Name of output file to write publication-ready PyMOL session file to. Should end in '.pse'; defaults to 'out.pse'.",
    required=False,
    default="out.pse",
)
@click.option(
    "-oi",
    "--outfile-interactive",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Name of output file to write ready-to-share interactive HTML file to. Should end in '.html'; defaults to 'out.html'.",
    required=False,
    default="out.html",
)
@click.option(
    "-fc",
    "--fitness-column",
    type=click.STRING,
    help="Name of the column in the fitness-file (-f/--fitness-file) that contains fitness values (e.g. LogEffect). If not defined, will default to 'fitness'.",
    required=False,
    default="fitness",
)
@click.option(
    "-ri",
    "--residue-index-column",
    type=click.STRING,
    help="Name of the column in the fitness-file (-f/--fitness-file) that contains residue indices (e.g. 1, 2, .. n). If not defined, will default to 'residue_index'.",
    required=False,
    default="residue_index",
)
@click.option(
    "-wt",
    "--wildtype-column",
    type=click.STRING,
    help="Name of the column in the fitness-file (-f/--fitness-file) that contains wildtype residues (e.g. L, G, N). If not defined, will default to 'wildtype'.",
    required=False,
    default="wildtype",
)
@click.option(
    "-mu",
    "--mutant-column",
    type=click.STRING,
    help="Name of the column in the fitness-file (-f/--fitness-file) that contains mutant residues (e.g. L, G, N). If not defined, will default to 'mutant'.",
    required=False,
    default="mutant",
)
@click.option(
    "-c",
    "--confidence-column",
    type=click.STRING,
    help="Name of the column in the fitness-file (-f/--fitness-file) that contains confidence values (e.g. counts). If not defined then LogoPlots in the HTML view will not display confidences.",
    required=False,
)
def render(
    pdb_file: Optional[str] = None,
    fitness_file: Optional[str] = None,
    fitness_threshold: Optional[float] = None,
    outfile_publication: Optional[str] = None,
    outfile_interactive: Optional[str] = None,
    fitness_column: Optional[str] = None,
    residue_index_column: Optional[str] = None,
    wildtype_column: Optional[str] = None,
    mutant_column: Optional[str] = None,
    confidence_column: Optional[str] = None,
):

    if not outfile_publication[-4:] == ".pse":
        raise ValueError("--op/--outfile-publication should end in '.pse'.")
    elif not outfile_interactive[-5:] == ".html":
        raise ValueError("--oi/--outfile-interactive should end in '.html'.")

    fitness_dict = FitnessFactory(
        fitness_file,
        resindex_colname=residue_index_column,
        wildtype_colname=wildtype_column,
        mutant_colname=mutant_column,
        fitness_colname=fitness_column,
        confidence_colname=confidence_column,
    ).get_fitness_basedict()
    complex = ComplexFactory(pdb_file).load_pdb()
    complex_rdkit = ComplexFactory(pdb_file).load_pdb_rdkit()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    PublicationView(
        filled_aligned_fitness_dict,
        complex,
        complex_rdkit,
        fitness_threshold=fitness_threshold,
        output_session_file=outfile_publication,
    ).render()

    # TODO: add more logging in the below Class.
    InteractiveView(
        filled_aligned_fitness_dict,
        complex,
        complex_rdkit,
        fitness_threshold=fitness_threshold,
        output_session_file=outfile_interactive,
    ).render()
