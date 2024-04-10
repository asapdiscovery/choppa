import MDAnalysis
from MDAnalysis.lib.util import NamedStream
from io import StringIO
import warnings

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