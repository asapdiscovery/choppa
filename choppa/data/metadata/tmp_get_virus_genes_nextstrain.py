from choppa.nextstrain import fetch_nextstrain_json
from choppa.data.metadata.resources import NEXTSTRAIN_METADATA
from urllib.parse import urlparse

import json


def get_url(virus):
    # same as in choppa.nextstrain but skips virus/gene name checks

    # get the URL from metadata for this virus
    parsed_url = urlparse(nextstrain_metadata[virus]["URL"])
    assert parsed_url.netloc == "nextstrain.org"

    # Extract the path
    path = parsed_url.path[1:]

    # return the properly formatted path to the mutation data on the database backend
    data_path = "_".join(path.split("/"))
    return f"https://data.nextstrain.org/{data_path}.json"


with open(NEXTSTRAIN_METADATA) as json_data:
    nextstrain_metadata = json.load(json_data)

for virus, data in nextstrain_metadata.items():
    list_of_genes = []
    for gene in list(
        fetch_nextstrain_json(get_url(virus))["meta"]["genome_annotations"].keys()
    ):
        if not gene == "nuc":
            list_of_genes.append(gene)
    print(virus, list_of_genes)
