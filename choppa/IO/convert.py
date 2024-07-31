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
    fitness_df = pd.DataFrame(
        json.load(open(json_file))["ZIKV NS2B-NS3 (Open)"]["mut_metric_df"]
    )

    if gene:
        print(f"Available genes: {set(fitness_df['gene'].values)}")

        fitness_df_for_gene = fitness_df[fitness_df["gene"] == gene]
        fitness_df_for_gene = fitness_df_for_gene[
            fitness_df_for_gene["site"].between(207, 379)
        ]
        if len(fitness_df_for_gene) == 0:
            raise ValueError(f"Gene '{gene}' not found in data:\n{fitness_df}")
        else:
            logger.info(
                f"Extracted {len(fitness_df_for_gene)} entries across {len(fitness_df_for_gene.groupby(by='site'))} sites during JSON->DataFrame conversion"
            )
            return fitness_df_for_gene
    else:
        return fitness_df


def ns2b3_reset_residcs(df):
    """ns2b3 has the same indices between the two chains. Super annoying, resetting that here."""
    new_idcs_col = []
    for idx in df["reference_site"].values:
        if "(NS2B) " in idx:
            new_idcs_col.append(idx.replace("(NS2B) ", ""))
        elif "(NS3) " in idx:
            new_idcs_col.append(130 + int(idx.replace("(NS3) ", "")))
    df["reference_site"] = new_idcs_col

    return df


if __name__ == "__main__":
    fitness_df = phylo_json_to_df(sys.argv[1])

    fitness_df = ns2b3_reset_residcs(fitness_df)

    fitness_df.to_csv(sys.argv[2], index=False)
