import json
import pandas as pd
import logging, sys

from choppa.data.toy_data.resources import TOY_PHYLO_DATA

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()


def phylo_json_to_df(json_file, gene=None):
    """
    Converts a dataset of phylogenetics fitness data formatted in JSON into
    a table format as DataFrame as used in `choppa.IO.input`. This is prone to
    breaking (depending on how the JSON is formatted). Ideally the JSON
    has been generated with https://github.com/jbloomlab/SARS2-mut-fitness/blob/main/scripts/export_fitness_to_json.py

    `gene` can be specified to only export a specific gene into the dataframe.
    """
    fitness_df = pd.DataFrame(json.load(open(json_file))["data"])

    if gene:
        fitness_df_for_gene = fitness_df[fitness_df["gene"] == gene]
        if len(fitness_df_for_gene) == 0:
            raise ValueError(f"Gene '{gene}' not found in data:\n{fitness_df}")
        else:
            logger.info(
                f"Extracted {len(fitness_df_for_gene)} entries across {len(fitness_df_for_gene.groupby(by='site'))} sites during JSON->DataFrame conversion"
            )
            return fitness_df_for_gene
    else:
        return fitness_df


def nextstrain_to_csv(nextstrain_tsv):
    """ """


if __name__ == "__main__":
    phylo_json_to_df(TOY_PHYLO_DATA, "nsp9")
