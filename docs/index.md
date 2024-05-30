
Welcome to choppa's documentation!
=========================================================

`Choppa` is a package for visualizing fitness data (e.g. deep mutational scanning or phylogenetics data) on 3D protein structures. It has been built as an end-to-end CLI tool that ingests a `PDB` file (protein structure) and a `CSV` file (fitness data) and produces two files:
- A `HTML` file with interactive features such as `logoplot` pop-ups when residues are hovered over.
- A `PyMOL` session file with the same coloring as the HTML file. This view can be ray-traced with PyMOL and used for publication-style figures.  

```{toctree}
:maxdepth: 5
:caption: Contents:
:hidden:

getting_started
API/index
```
