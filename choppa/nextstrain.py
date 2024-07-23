import re
import math
import copy
import click
import requests
import pandas as pd
import altair as alt
from Bio import Phylo
from urllib.parse import urlparse
from augur.utils import annotate_parents_for_tree


"""
CLI endpoint to scrape nextstrain, take in virus name and gene and whether csv should be written to file
in metadata have dict of viruses and URLS and genes per virus
take Will's steps
return DF and write to file if requested
"""

