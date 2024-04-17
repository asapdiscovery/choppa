import MDAnalysis
from MDAnalysis.lib.util import NamedStream
from MDAnalysis.analysis import distances
from MDAnalysis.core.groups import AtomGroup
import pymol2
from Bio.PDB.PDBIO import PDBIO

from io import StringIO
import warnings
import tempfile

def get_ligand_resnames_from_pdb_str(PDB_str, remove_solvent=True):
    """
    Uses MDAnalysis to figure out what residue names the ligand(s) in the protein PDB (str) has/have.

    Uses StringIO to circumvent having to write to memory.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # hides MDA RunTimeWarning that complains about string IO
        u = MDAnalysis.Universe(NamedStream(StringIO(PDB_str), "complex.pdb"))

    if remove_solvent:
        ag = u.select_atoms("not protein and not (name H* or type OW)")
    else:
        ag = u.select_atoms("not protein")
    resnames = set(ag.resnames)
    return list(resnames)

def biopython_to_mda(BP_complex):
    """
    Converts a biopython protein object to an MDAnalysis one.
    """
    with tempfile.TemporaryDirectory() as tmpdirname:
        io=PDBIO()
        io.set_structure(BP_complex)
        io.save(f"{tmpdirname}/tmp_while_hmo_helps_write_to_stream.pdb")

        u = MDAnalysis.Universe(f"{tmpdirname}/tmp_while_hmo_helps_write_to_stream.pdb")
    return u   

def get_pdb_components(PDB_str, remove_solvent=True):
    """
    Split a protein-ligand pdb into protein and ligand components
    :param PDB_str:
    :return:
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # hides MDA RunTimeWarning that complains about string IO
        u = MDAnalysis.Universe(NamedStream(StringIO(PDB_str), "complex.pdb"))

    if remove_solvent:
        ag = u.select_atoms("not (name H* or type OW)")
    
    ligand = u.select_atoms("not protein")
    protein = u.select_atoms("protein")

    return ligand, protein


def process_ligand(ligand):
    """
    Add bond orders to a pdb ligand in an MDA universe object.
    1. load PDB into PyMol session (PyMOL does the bond guessing)
    2. write ligand to stream as SDF
    3. Read the stream into an RDKit molecule
    """
    with tempfile.TemporaryDirectory() as tmpdirname:
        ligand.write(f"{tmpdirname}/lig_tmp_while_hmo_helps_write_to_stream.pdb")

        p = pymol2.PyMOL()
        p.start()
        p.cmd.load(f"{tmpdirname}/lig_tmp_while_hmo_helps_write_to_stream.pdb")
        p.cmd.save(f"{tmpdirname}/lig_tmp_while_hmo_helps_write_to_stream.sdf", "all", 0) # writes all states, so should be able to handle multi-ligand
        p.stop()

        with open(f"{tmpdirname}/lig_tmp_while_hmo_helps_write_to_stream.sdf","r") as f:
            string = f.read()
    return string

def process_protein(protein):
    """
    Returns the string for the protein in an MDA universe object.
    """
    with tempfile.TemporaryDirectory() as tmpdirname:
        protein.write(f"{tmpdirname}/prot_tmp_while_hmo_helps_write_to_stream.pdb")
        with open(f"{tmpdirname}/prot_tmp_while_hmo_helps_write_to_stream.pdb","r") as f:
            string = f.read()
    return string

def split_pdb_str(PDB_str):
    """
    From a PDB string, gets the string for the protein and (if present) the ligand SDF (with guessed
    bond orders).

    Inspired by https://gist.github.com/PatWalters/c046fee2760e6894ed13e19b8c99193b
    TODO: set below functions through NamedStream instead of tmpdir
    """
    ligand_pdb, protein_pdb = get_pdb_components(PDB_str)
    return process_ligand(ligand_pdb),  process_protein(protein_pdb)

def show_contacts(
    pymol_instance,
    selection_residues,
    selection_lig,
    contact_color,
    bigcutoff=4.0,
):
    """
    Heavily reduced PyMOL plugin that provides show_contacts command and GUI for highlighting good and bad 
    polar contacts. Factored out of clustermols by Matthew Baumgartner. 

    Returns:
    List of contacts
    """
    # if the group of contacts already exist, delete them
    pymol_instance.cmd.delete(pymol_instance.cmd.get_legal_name("contacts"))

    # ensure only N and O atoms are in the selection
    all_don_acc_lig = selection_lig + " and (donor or acceptor)"
    all_don_acc_res = selection_residues + " and  (donor or acceptor) and sidechain"
    all_don_acc_res_backbone =  selection_residues + " and  (donor or acceptor) and backbone"

    # if theses selections turn out not to have any atoms in them, pymol throws cryptic errors when calling the dist function like:
    # 'Selector-Error: Invalid selection name'
    all_lig_sele_count = pymol_instance.cmd.select("all_don_acc_lig_sele", all_don_acc_lig)
    all_res_sele_count = pymol_instance.cmd.select("all_don_acc_res_sele", all_don_acc_res)
    _ = pymol_instance.cmd.select("all_don_acc_res_bb_sele", all_don_acc_res_backbone) # not sure if we need to check for this backbone one?
    if not all_lig_sele_count and all_res_sele_count:
        return False
    
    # now make the actual contacts
    # first with sidechains, color as requested
    contacts_name = f"{selection_residues}_contacts_sc"
    pymol_instance.cmd.distance(
        contacts_name, "all_don_acc_lig_sele", "all_don_acc_res_sele", bigcutoff, mode=0
    )
    pymol_instance.cmd.set("dash_radius", "0.09", contacts_name)
    pymol_instance.cmd.set("dash_color", contact_color, contacts_name)
    pymol_instance.cmd.hide("labels", contacts_name)
    
    # now contacts with backbone, color these green
    contacts_name = f"{selection_residues}_contacts_bb"
    pymol_instance.cmd.distance(
        contacts_name, "all_don_acc_lig_sele", "all_don_acc_res_bb_sele", bigcutoff, mode=0
    )
    pymol_instance.cmd.set("dash_radius", "0.09", contacts_name)
    pymol_instance.cmd.set("dash_color", "green", contacts_name)
    pymol_instance.cmd.hide("labels", contacts_name)

    return True

def get_contacts_mda(
    complex,
    bigcutoff=4.1, # for some reason MDA needs a little 0.1A nudge to agree with PyMOL measurements 
    remove_solvent=True,
):
    """
    Use MDAnalysis to generate a dictionary of distance endpoint xyz coordinates between atoms in the ligand
    and protein residues.
    """
    contacts = []
    u = biopython_to_mda(complex)
    if remove_solvent:
        u = u.select_atoms("not (name H* or type OW)")
    
    lig = u.select_atoms("not protein")
    prot = u.select_atoms("protein")

    distances_array = distances.distance_array(AtomGroup(lig), AtomGroup(prot))
    distances_array_within_cutoff = distances_array <= bigcutoff # makes the NxM array be boolean based on cutoff

    for contacted_lig_at, protein_sequence_distance_bool in zip(lig, distances_array_within_cutoff):
        # for this atom in the ligand, find any protein atoms that are True (i.e. below cutoff)
        contacted_res_ats = [prot[i] for i, contact in enumerate(protein_sequence_distance_bool) if contact]

        for contacted_res_at in contacted_res_ats:
            # we only want to show HBD/HBA
            if contacted_lig_at.type == "N" and contacted_res_at.type == "O" or \
                            contacted_lig_at.type == "O" and contacted_res_at.type == "N":
                contacts.append([contacted_lig_at, contacted_res_at])
    
    return contacts
        

        
            