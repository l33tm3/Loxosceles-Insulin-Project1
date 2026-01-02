
import os
import math
import numpy as np
from Bio import PDB
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser

def get_residue_center(residue):
    """Calculate Geometric Center of a residue."""
    atoms = [a.get_vector().get_array() for a in residue if a.element != 'H']
    return np.mean(atoms, axis=0) if atoms else None

def measure_distance():
    print("--- Bioinformatic Distance Measurement Tool ---")
    
    # Files
    toxin_file = r"c:\Users\nanos\OneDrive\Documentos\Curitas\fold_rinconsito_1_model_0.cif"
    pdb_1irk = "pdb1irk.ent" # Standard download name
    
    # 1. Load Toxin
    if not os.path.exists(toxin_file):
        print(f"Error: {toxin_file} not found.")
        return
    
    parser_cif = MMCIFParser(QUIET=True)
    toxin_struct = parser_cif.get_structure("Toxin", toxin_file)
    
    # 2. Identify Toxin Active Site (Based on analysis: HIS 38 + GLU 58)
    # Note: Toxin numbering is model-specific.
    toxin_res_targets = [38, 58] 
    toxin_centers = []
    
    print(f"\n[Toxin] Target Residues (Potential Active Site): {toxin_res_targets}")
    for model in toxin_struct:
        for chain in model:
            for res in chain:
                if res.get_id()[1] in toxin_res_targets:
                    center = get_residue_center(res)
                    toxin_centers.append((res.get_resname(), res.get_id()[1], center))
                    print(f"  Found {res.get_resname()} {res.get_id()[1]} at {center}")
                    
    if not toxin_centers:
        print("  Warning: Targeted toxin residues not found.")

    # 3. Load 1IRK
    start_dir = os.getcwd()
    pdbl = PDB.PDBList()
    # Check if 1IRK is present or download
    try:
        if not os.path.exists(pdb_1irk):
             print("\n[Downloading] 1IRK from PDB...")
             pdbl.retrieve_pdb_file("1IRK", pdir=".", file_format="pdb")
    except Exception as e:
        print(f"Warning: Could not download 1IRK: {e}")

    parser_pdb = PDBParser(QUIET=True)
    try:
        irk_struct = parser_pdb.get_structure("1IRK", pdb_1irk)
        print("\n[1IRK] Loaded Structure.")
    except:
        print("  Error: Could not load 1IRK structure.")
        return

    # 4. Identify 1IRK "Sugar Receptor" / Active Site
    # 1IRK is a Kinase. Active site: Asp 1132 (Catalytic Loop), Lys 1030
    irk_targets = [1132, 1030] 
    irk_centers = []
    
    print(f"[1IRK] Target Residues (Kinase Active Site): {irk_targets}")
    for model in irk_struct:
        for chain in model:
            for res in chain:
                if res.get_id()[1] in irk_targets:
                    center = get_residue_center(res)
                    irk_centers.append((res.get_resname(), res.get_id()[1], center))
                    print(f"  Found {res.get_resname()} {res.get_id()[1]} at {center}")

    # 5. Measure Distances
    print("\n--- Distances (Angstroms) ---")
    print("Note: Unless structures are superimposed, these distances are relative to their file origins.")
    
    if toxin_centers and irk_centers:
        for t_name, t_id, t_vec in toxin_centers:
            for i_name, i_id, i_vec in irk_centers:
                dist = np.linalg.norm(t_vec - i_vec)
                print(f"Distance {t_name}{t_id} (Toxin) <--> {i_name}{i_id} (1IRK): {dist:.2f} A")
    else:
        print("Could not calculate distances (missing residues).")

if __name__ == "__main__":
    measure_distance()
