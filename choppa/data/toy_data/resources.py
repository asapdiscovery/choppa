import pkg_resources

TOY_COMPLEX = pkg_resources.resource_filename(
    __name__, "toy_complex_SARS-CoV-2-Mac1.pdb"
)

TOY_FITNESS_DATA_COMPLETE = pkg_resources.resource_filename(
    __name__, "toy_fitness_data_01_complete.csv"
)

TOY_FITNESS_DATA_COMPLETE_NOCONF = pkg_resources.resource_filename(
    __name__, "toy_fitness_data_01_complete_noconf.csv"
)

TOY_FITNESS_DATA_TRUNCATED = pkg_resources.resource_filename(
    __name__, "toy_fitness_data_02_truncated_25.csv"
)

TOY_FITNESS_DATA_SECTIONED = pkg_resources.resource_filename(
    __name__, "toy_fitness_data_03_sectioned.csv"
)
TOY_HYLO_DATA = pkg_resources.resource_filename(
    __name__, "toy_phylo_fitness_data.json"
)
