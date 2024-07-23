import re
import math
import copy
import requests
import json
from urllib.parse import urlparse

import pandas as pd
import click

import altair as alt
from Bio import Phylo
from augur.utils import annotate_parents_for_tree

from choppa.data.metadata.resources import NEXTSTRAIN_METADATA

AMINO_ACIDS = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "*",
]


def validate_virus_gene(virus, gene):
    """Checks that the provided virus and gene strings are available in NextStrain."""
    with open(NEXTSTRAIN_METADATA) as json_data:
        nextstrain_metadata = json.load(json_data)

    if not virus in nextstrain_metadata.keys():
        raise ValueError(
            f"Virus named '{virus}' not in viruses available on NextStrain: {list(nextstrain_metadata.keys())}"
        )

    if not gene in nextstrain_metadata[virus]["genes"]:
        raise ValueError(
            f"Gene named '{gene}' not in genes on virus '{virus}' available on NextStrain: {nextstrain_metadata[virus]['genes']}"
        )
    return nextstrain_metadata


def get_url(virus, gene):
    """
    Returns the full download URL for downstream GET requests.
    """
    # validate first
    nextstrain_metadata = validate_virus_gene(virus, gene)

    # get the URL from metadata for this virus
    parsed_url = urlparse(nextstrain_metadata[virus]["URL"])
    assert parsed_url.netloc == "nextstrain.org"

    # Extract the path
    path = parsed_url.path[1:]

    # return the properly formatted path to the mutation data on the database backend
    data_path = "_".join(path.split("/"))
    return (
        f"https://data.nextstrain.org/{data_path}.json",
        nextstrain_metadata[virus]["URL"],
    )


def fetch_nextstrain_json(url):
    # Make a GET request to the specified URL
    response = requests.get(url)

    # Raise an exception if the request failed
    response.raise_for_status()

    # Parse the JSON data and return it
    return response.json()


def fetch_nextstrain_root_sequence(url):
    # Header to request the root sequence data
    headers = {"Accept": "application/vnd.nextstrain.dataset.root-sequence+json"}

    # Make a GET request to the specified URL
    response = requests.get(url, headers=headers)

    # Raise an exception if the request failed
    response.raise_for_status()
    # raise specific errors here?

    return response.json()


def nextstrain_json_to_tree(json_dict, root=True, parent_cumulative_branch_length=None):
    """
    Follow the exact same logic that the NextStrain backend follows when calculating per-residue mutation frequencies.
    By John Huddleston.
    """
    # Check for v2 JSON which has combined metadata and tree data.
    if root and "meta" in json_dict and "tree" in json_dict:
        json_dict = json_dict["tree"]

    node = Phylo.Newick.Clade()

    # v1 and v2 JSONs use different keys for strain names.
    if "name" in json_dict:
        node.name = json_dict["name"]
    else:
        node.name = json_dict["strain"]

    # Assign all non-children attributes.
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    # Only v1 JSONs support a single `attr` attribute.
    if hasattr(node, "attr"):
        node.numdate = node.attr.get("num_date")
        node.cumulative_branch_length = node.attr.get("div")

        if "translations" in node.attr:
            node.translations = node.attr["translations"]
    elif hasattr(node, "node_attrs"):
        node.cumulative_branch_length = node.node_attrs.get("div")

    node.branch_length = 0.0
    if parent_cumulative_branch_length is not None and hasattr(
        node, "cumulative_branch_length"
    ):
        node.branch_length = (
            node.cumulative_branch_length - parent_cumulative_branch_length
        )

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [
            nextstrain_json_to_tree(
                child,
                root=False,
                parent_cumulative_branch_length=node.cumulative_branch_length,
            )
            for child in json_dict["children"]
        ]

    if root:
        node = annotate_parents_for_tree(node)

    return node


def extract_tree_data(tree, attributes=None, include_internal_nodes=False):
    """
    Further process mutation data into a usable dataframe.
    By John Huddleston.
    """
    records = []

    if attributes is None:
        attributes = sorted(
            set(tree.root.node_attrs.keys()) | set(tree.root.branch_attrs.keys())
        )

    for node in tree.find_clades():
        if node.is_terminal() or include_internal_nodes:
            record = {"name": node.name}

            for attribute in attributes:
                if attribute in node.node_attrs:
                    value = node.node_attrs[attribute]
                elif attribute in node.branch_attrs:
                    value = node.branch_attrs[attribute]
                else:
                    print(
                        f"Could not find attribute '{attribute}' for node '{node.name}'."
                    )
                    value = None

                if value is not None:
                    if isinstance(value, dict) and "value" in value:
                        value = value["value"]

                record[attribute] = value

            records.append(record)

    # Convert records to a data frame
    return pd.DataFrame(records)


def parse_mutations(mutation):
    """Regex for parsing mutations into (from, position, to)"""
    pattern = re.compile(r"([^0-9]+)(\d+)([^0-9]+)")
    match = pattern.search(mutation)
    if match:
        return match.groups()
    else:
        return None


# def calculate_entropy(tree, mutations, gene):
# This is not needed for now, keeping here for potential later use
#     # The mutation at the root of the tree (ancestor) for each position
#     ancestral_state = dict()
#     # Observed count of aa's for each position
#     counts = dict()
#     # Count of total tips
#     tips = 0

#     # Recursive function to traverse the tree
#     def recurse_tree(node, state):
#         nonlocal tips
#         nonlocal ancestral_state
#         nonlocal counts

#         # Retrieve and parse the mutations for this node
#         mutations_at_node = mutations.get(node.name, {})
#         # Update the state and the ancestral 'sequence'
#         for cds, mutation_list in mutations_at_node.items():
#             if cds == gene:
#                 for mutation in mutation_list:
#                     aa_from, position, aa_to = parse_mutations(mutation)
#                     if position not in ancestral_state.keys():
#                         ancestral_state[position] = aa_from
#                     state[position] = aa_to

#         # If the node doesn't have children, update the count
#         if not node.clades:
#             tips += 1
#             for position, aa in state.items():
#                 if position not in counts:
#                     counts[position] = {}
#                     counts[position][aa] = 1
#                 elif aa not in counts[position]:
#                     counts[position][aa] = 1
#                 else:
#                     counts[position][aa] += 1
#         else:
#             for clade in node.clades:
#                 recurse_tree(clade, copy.deepcopy(state))

#     # Start recursion at the root
#     recurse_tree(tree.root, {})

#     assert tips == len(tree.get_terminals())

#     # Calculate the entropy for each position
#     entropy = []
#     m = 0
#     i = 0
#     for position in counts.keys():
#         nObserved = 0
#         for observedCount in list(counts[position].values()):
#             nObserved += observedCount
#         nUnobserved = tips - nObserved
#         if nUnobserved > 0:
#             if ancestral_state[position] in counts[position].keys():
#                 counts[position][ancestral_state[position]] += nUnobserved
#             else:
#                 counts[position][ancestral_state[position]] = nUnobserved

#         s = 0
#         for count in list(counts[position].values()):
#             a = count / tips
#             s += -1 * a * math.log(a)
#         if s > m:
#             m = s
#         entropy.append(
#             {"entropy": round(s, 3), "codon_position": int(position), "gene": gene}
#         )
#         i += 1

#     return pd.DataFrame(entropy)

# code for plotting an interactive sequence-mutation HTML figure like on NextStrain web UI:
# chart = (
#     alt.Chart(entropy_df)
#     .mark_bar()
#     .encode(
#         x="codon_position:Q",
#         y="entropy:Q",
#         tooltip=["codon_position:Q", "entropy:Q", "gene"],
#         color=alt.value("skyblue"),
#     )
#     .properties(width=1600, height=200, title=f"Entropy of {gene}")
# )
# brush = alt.selection_interval(encodings=["x"])  # Brush for selection
# overview = (
#     alt.Chart(entropy_df)
#     .mark_rect()
#     .encode(x="codon_position:Q", color=alt.value("skyblue"))
#     .mark_rect()
#     .add_params(brush)
#     .properties(width=1600, height=25, title="site zoom bar")
# )
# chart = chart.transform_filter(brush)
# combined_chart = alt.vconcat(chart, overview)
# combined_chart.save(f"{data_path}_{gene}_entropy.html")


def count_mutations_events(metadata_df, gene):
    mutations = pd.Series(
        metadata_df.mutations.values, index=metadata_df.name
    ).to_dict()
    # Count the number of times a mutation occurs on the tree
    mutation_counts = dict()
    for mutation_list in mutations.values():
        if gene in mutation_list.keys():
            for mutation in mutation_list[gene]:
                _, position, aa_to = parse_mutations(mutation)
                if position not in mutation_counts:
                    mutation_counts[position] = {}
                    mutation_counts[position][aa_to] = 1
                elif aa_to not in mutation_counts[position]:
                    mutation_counts[position][aa_to] = 1
                else:
                    mutation_counts[position][aa_to] += 1

    # Convert the counts to a data frame
    rows = []
    for position, mutations in mutation_counts.items():
        for mutation, count in mutations.items():
            rows.append((position, mutation, count))
    mutation_count_df = pd.DataFrame(rows, columns=["position", "mutation", "count"])
    mutation_count_df["position"] = mutation_count_df["position"].astype(int)
    mutation_count_df = mutation_count_df.sort_values(by="position").reset_index(
        drop=True
    )

    return mutation_count_df


def finalize_dataframe(mutation_count_df, root_sequence_json, outfile):
    root_sequence_df = pd.DataFrame(
        [(i + 1, aa) for i, aa in enumerate(root_sequence_json["ORF7a"])],
        columns=["position", "residue"],
    )

    # add them together so that we have a df with wildtypes, mutations and each mutation's count
    counts_df = root_sequence_df.merge(mutation_count_df, on="position")

    # now construct a choppa-style df. We need to iterate all possible mutations and note the count for each residue
    choppa_nextstrain_data = []
    for resi in range(1, len(root_sequence_json["ORF7a"])):
        # get the mutations for this residue number. This dataframe may be empty if there are no mutations in NextStrain.
        recorded_mutations = counts_df[counts_df["position"] == resi]

        if len(recorded_mutations) == 0:
            # easy: this residue has no mutants in NextStrain so just set frequencies to 0.
            for aa in AMINO_ACIDS:
                choppa_nextstrain_data.append(
                    {
                        "residue_index": resi,
                        "wildtype": root_sequence_json["ORF7a"][resi],
                        "mutant": aa,
                        "frequency": 0,
                    }
                )
        else:
            # we need to include the mutation frequencies found in NextStrain, then for the remaining
            # mutations we need to set frequencies to 0.
            for aa in AMINO_ACIDS:
                if aa not in recorded_mutations["mutation"].values:
                    frequency = 0  # this possible mutation isn't recorded in NextStrain
                else:
                    frequency = recorded_mutations[
                        recorded_mutations["mutation"] == aa
                    ]["count"].values[
                        0
                    ]  # but this one is

                choppa_nextstrain_data.append(
                    {
                        "residue_index": resi,
                        "wildtype": root_sequence_json["ORF7a"][resi],
                        "mutant": aa,
                        "frequency": frequency,
                    }
                )
    # add all rows into a single dataframe ready for usage with choppa.
    choppa_nextstrain_df = pd.DataFrame(choppa_nextstrain_data)
    if outfile:
        choppa_nextstrain_df.to_csv(outfile)
    return choppa_nextstrain_df
