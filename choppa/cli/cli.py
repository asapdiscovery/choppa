import click
import shutil
from typing import Optional

from choppa.IO.input import FitnessFactory, ComplexFactory
from choppa.align import AlignFactory
from choppa.render import PublicationView, InteractiveView


from choppa.data.metadata.resources import (
    MERS_REAL_GENE_TO_COMMUNITY_GENE,
    MERS_TARGET_TO_REAL_GENE,
)

from choppa.cli.utils import SpecialHelpOrder
from pathlib import Path


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
@click.option(
    "-pa",
    "--color-per-atom",
    is_flag=True,
    help="Color the fitness protein surface by atom rather than by residue: backbone atoms will always have a white fitness surface color.",
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
    color_per_atom: Optional[bool] = False,
):

    # check extensions
    if not Path(outfile_publication).suffix == ".pse":
        raise ValueError("--op/--outfile-publication should end in '.pse'.")

    if not Path(outfile_interactive).suffix == ".html":
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
        color_per_atom=color_per_atom,
        output_session_file=outfile_interactive,
    ).render()


@cli.command(
    name="nextstrain",
    help=". ",
    short_help="From the database of NextStrain-maintained pathogen analyses (https://nextstrain.org), generate a data format suitable for choppa.render.",
)
@click.option(
    "-v",
    "--virus",
    type=click.STRING,
    help="Name of the virus to download mutation data for. See https://nextstrain.org/pathogens for a list of available viruses.",
    required=True,
)
@click.option(
    "-g",
    "--gene",
    type=click.STRING,
    help="Name of the gene to download mutation data for. See e.g. https://nextstrain.org/zika for a view of available genes.",
    required=True,
)
@click.option(
    "-o",
    "--outfile",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Name of output file to write mutation data to in CSV. Should end in '.csv'.",
    required=True,
)
def nextstrain(
    virus: Optional[str] = None,
    gene: Optional[str] = None,
    outfile: Optional[str] = None,
):
    from choppa.nextstrain import (
        get_url,
        fetch_nextstrain_json,
        fetch_nextstrain_root_sequence,
        fetch_nextstrain_json_mers_cov,
        nextstrain_json_to_tree,
        extract_tree_data,
        count_mutations_events,
        finalize_dataframe,
    )

    # check extension
    if not Path(outfile).suffix == ".csv":
        raise ValueError("-o/--outfile should end in '.csv'.")

    download_url, nextstrain_tree_url = get_url(virus, gene)

    # Fetch the JSON data from the data URL
    # Fetch the root sequence data
    if virus == "MERS-CoV":
        tree_json = fetch_nextstrain_json_mers_cov()
        root_sequence_json = fetch_nextstrain_root_sequence(
            nextstrain_tree_url, MERS_COV=True
        )
    else:
        tree_json = fetch_nextstrain_json(download_url)
        root_sequence_json = fetch_nextstrain_root_sequence(nextstrain_tree_url)

    if root_sequence_json is None:
        # Fallback to tree_json if the root sequence is not available via the URL
        if "root_sequence" in tree_json:
            root_sequence_json = tree_json["root_sequence"]
        else:
            # Fail if no root sequence is available
            raise ValueError(
                "Root sequence is missing from the Nextstrain API and the main tree data."
            )

    # Make a tree from the JSON data
    tree = nextstrain_json_to_tree(tree_json)

    # Extract the mutations from the tree
    metadata_df = extract_tree_data(
        tree, attributes=["mutations"], include_internal_nodes=True
    )

    # Count terminal mutations
    treating_mers = False
    if virus == "MERS-CoV":
        # rename to match inconsistent gene naming in this community contribution NextStrain
        gene = MERS_REAL_GENE_TO_COMMUNITY_GENE[MERS_TARGET_TO_REAL_GENE[gene]]
        treating_mers = True

    mutation_count_df = count_mutations_events(metadata_df, gene)

    # finalize dataframe by adding mutations and root sequence together
    _ = finalize_dataframe(
        mutation_count_df, root_sequence_json, gene, outfile, mers=treating_mers
    )
