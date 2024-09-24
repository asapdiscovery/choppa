import pkg_resources

# Fitness data in JSON format
NEXTSTRAIN_METADATA = pkg_resources.resource_filename(__name__, "nextstrain.json")

# MERS is a community dataset, so not directly queriable. Providing JSON here, download
# from https://raw.githubusercontent.com/wqshi/merse/refs/heads/master/auspice/merse_tree.json
MERS_GENOME_TREE_JSON = pkg_resources.resource_filename(
    __name__, "MERS-CoV_auspice_tree.json"
)

# MERS is a community dataset and the gene names are borked. We'll rename them here manually.
# curated manually using FEATURES.CDS indices in MERS-CoV-reference.gb
MERS_REAL_GENE_TO_COMMUNITY_GENE = {
    "S": "G128_gp02",
    "orf3": "G128_gp03",
    "orf4a": "G128_gp04",
    "orf4b": "G128_gp05",
    "orf5": "G128_gp06",
    "E": "G128_gp07",
    "M": "G128_gp08",
    "N": "G128_gp09",
    "orf8b": "G128_gp10",
}

MERS_TARGET_TO_REAL_GENE = {
    "S": "S",
    "ORF3": "orf3",
    "ORF4a": "orf4a",
    "ORF4b": "orf4b",
    "ORF5": "orf5",
    "E": "E",
    "M": "M",
    "N": "N",
    "ORF8b": "orf8b",
}

# MERS root sequence - not provided in community JSON so had to make one
# ourselves based on https://raw.githubusercontent.com/nextstrain/mers/refs/heads/master/config/mers_reference.gb
MERS_GENOME_ROOT_SEQUENCE = pkg_resources.resource_filename(
    __name__, "MERS-CoV-reference.gb"
)
