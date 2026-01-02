
import os
import numpy as np
from Bio import PDB
from Bio.PDB import Superimposer, MMCIFParser, PDBParser, PDBIO

def perform_superposition():
    print("--- Structural Superposition ---")
    
    # Paths
    toxin_file = r"c:\Users\nanos\OneDrive\Documentos\Curitas\fold_rinconsito_1_model_0.cif"
    pdb_1irk = "pdb1irk.ent"
    
    # Check Files
    if not os.path.exists(toxin_file) or not os.path.exists(pdb_1irk):
        print("Error: Input files missing.")
        return

    # Load Structures
    parser_cif = MMCIFParser(QUIET=True)
    toxin_structure = parser_cif.get_structure("Toxin", toxin_file)
    
    parser_pdb = PDBParser(QUIET=True)
    irk_structure = parser_pdb.get_structure("1IRK", pdb_1irk)

    # Convert to atom lists
    # We need to select the SPECIFIC atoms to align.
    # Toxin: HIS 38 (CA, CB..), GLU 58 (CA, CB..)
    # 1IRK: LYS 1030, ASP 1132
    
    # Helper to find atoms
    def get_atoms(structure, chain_id, res_id):
        # Note: simplistic lookup, assuming first model and chain if chain_id None
        for model in structure:
            for chain in model:
                if chain_id and chain.id != chain_id: continue
                for res in chain:
                    if res.get_id()[1] == res_id:
                        # Return CA atoms for backbone alignment or all common atoms
                        return [res['CA']] # Aligning Alpha Carbons for simplicity/robustness first
        return []

    # Select Atoms for Alignment
    # Toxin: Res 38, 58
    fixed_atoms = [] # 1IRK (Target)
    moving_atoms = [] # Toxin (Mobile)
    
    # 1IRK Targets (CA atoms)
    # Note: 1IRK might have multiple chains, usually chain A is the kinase domain
    atom_irk_1 = get_atoms(irk_structure, None, 1030) # LYS
    atom_irk_2 = get_atoms(irk_structure, None, 1132) # ASP
    
    # Toxin Targets (CA atoms)
    atom_tox_1 = get_atoms(toxin_structure, None, 38) # HIS
    atom_tox_2 = get_atoms(toxin_structure, None, 58) # GLU
    
    if atom_irk_1 and atom_irk_2 and atom_tox_1 and atom_tox_2:
        fixed_atoms = atom_irk_1 + atom_irk_2
        moving_atoms = atom_tox_1 + atom_tox_2
    else:
        print("Error: Could not locate all target atoms for alignment.")
        return

    # Superimposer
    super_imposer = Superimposer()
    super_imposer.set_atoms(fixed_atoms, moving_atoms)
    
    print(f"Alignment calculated on {len(fixed_atoms)} atom pairs (Active Site CA).")
    print(f"RMSD: {super_imposer.rms:.4f}")

    # Apply transformation to the whole Toxin structure
    super_imposer.apply(toxin_structure.get_atoms())
    
    # Save aligned toxin
    io = PDBIO()
    io.set_structure(toxin_structure)
    io.save("toxin_aligned_to_1irk.pdb")
    print("Saved aligned structure to: toxin_aligned_to_1irk.pdb")

if __name__ == "__main__":
    perform_superposition()
