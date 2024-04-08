from Bio.PDB import PDBParser
import csv
import random

TOY_FITNESS_THRESHOLD=0.9

def write_toy_fitness_csv(toy_fitness_data, out_name):
    """
    """
    with open(out_name, "w") as writefile:
        writer = csv.writer(writefile)
        writer.writerow(["residue_index", "wildtype", "mutant", "fitness", "confidence"])
        for row in toy_fitness_data:
            writer.writerow(row)

    return out_name

def create_toy_fitness_csv(add_confidence=True):
    one_letter ={
                'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',
                'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',
                'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',
                'GLY':'G', 'PRO':'P', 'CYS':'C'
                }

    # load the toy protein-ligand complex
    structure = PDBParser(QUIET=True).get_structure("Mac1", "toy_complex_SARS-CoV-2-Mac1.pdb")

    # iterate over the complex's residues
    toy_fitness_data = []
    for i, res in enumerate(structure.get_residues(), start=1):
        resname_3 = res.get_resname()
        if resname_3 not in one_letter.keys():
            continue
        resname_1 = one_letter[resname_3]
        
        # find out what residues this 'wildtype' residue could mutate into
        mutants = [ v for _,v in one_letter.items() if not v == resname_1]

        # also add the stop codon
        mutants += ["X"]

        # add wildtype to toy fitness data
        toy_fitness_data.append([
                i, # we could use the actual PDB indexing (`res.get_id()[1]`), but fitness data never does this so best to stick to [1..n]
                resname_1, # this is 'wildtype'
                resname_1,
                1.0, # fitness value
                int(random.uniform(1, 5000)) # fitness confidence, in real world would be something like `number of reads in well`
            ])
        
        # add mutants to toy fitness data
        for mut in mutants:
            toy_fitness_data.append([
                i, 
                resname_1, 
                mut,
                round(random.uniform(-5.0, 1.0), 2), 
                int(random.uniform(1, 5000)) 
            ])
        
    if not add_confidence:
        # just remove it.
        toy_fitness_data = [ dat[:4] for dat in toy_fitness_data ]

    return toy_fitness_data, structure


# write out the toy fitness dataset
toy_fitness_data, _ = create_toy_fitness_csv()
print("Full fitness set:", len(toy_fitness_data))
write_toy_fitness_csv(toy_fitness_data, "toy_fitness_data_01_complete.csv")

# write out the toy fitness dataset without confidence column
toy_fitness_data, _ = create_toy_fitness_csv(add_confidence=False)
print("Full fitness set without confidence:", len(toy_fitness_data))
write_toy_fitness_csv(toy_fitness_data, "toy_fitness_data_01_complete_noconf.csv")

# also make a toy fitness dataset that is only a section of the toy PDB so we can test alignment
toy_fitness_data, structure = create_toy_fitness_csv()
truncate_length = 25
start_index = truncate_length # we'll just shave off a few bits at the start and end.
end_index = len(list(structure.get_residues())) - truncate_length
truncated_toy_fitness_data = [ dat for dat in toy_fitness_data if start_index < dat[0] < end_index ]
print("Trunacted fitness set:", len(truncated_toy_fitness_data))
write_toy_fitness_csv(truncated_toy_fitness_data, f"toy_fitness_data_02_truncated_{truncate_length}.csv")

# also make a test case with multiple fragmented sections to further test alignment
toy_fitness_data, _ = create_toy_fitness_csv()
sections = list(range(5, 35)) + list(range(82, 110)) + list(range(130, 162))
sectioned_toy_fitness_data = [ dat for dat in toy_fitness_data if dat[0] in sections ]
print("Sectioned fitness set:", len(sectioned_toy_fitness_data))
write_toy_fitness_csv(sectioned_toy_fitness_data, f"toy_fitness_data_03_sectioned.csv")
