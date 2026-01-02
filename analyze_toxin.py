
import os
from Bio import PDB
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import Superimposer

def analyze_structure():
    # Paths
    toxin_file = r"c:\Users\nanos\OneDrive\Documentos\Curitas\fold_rinconsito_1_model_0.cif"
    
    if not os.path.exists(toxin_file):
        print(f"Error: File not found at {toxin_file}")
        return

    # Parse Toxin Structure
    parser = MMCIFParser()
    structure = parser.get_structure("Toxin", toxin_file)
    
    print("Analyzing Toxin Structure...")
    catalytic_candidates = []
    
    # Identify Catalytic Residues (His, Asp, Glu)
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in ["HIS", "ASP", "GLU"]:
                    catalytic_candidates.append(residue)
                    print(f"Found Candidate: {residue.get_resname()} {residue.get_id()[1]}")

    # Download 1IRK
    pdbl = PDB.PDBList()
    try:
        # Check if file exists to avoid redownloading if run multiple times
        # PDBList downloads to current dir by default in lowercase
        pdbl.retrieve_pdb_file("1IRK", pdir=".", file_format="pdb")
        pdb_file = "pdb1irk.ent" # Default filename format from retrieve_pdb_file
        
        if not os.path.exists(pdb_file):
             # Try alternative name if needed, or check where it downloaded
             print("Warning: Standard 1IRK file not found immediately, checking typical paths...")
             # (In a real scenario we might need to search, but let's assume standard behavior for now)

        pdb_parser = PDB.PDBParser()
        # Loading 1IRK
        structure_1irk = pdb_parser.get_structure("1IRK", pdb_file)
        print("\nSuccessfully loaded 1IRK.")

        # Minimal superimposition check (just to verify we can align)
        # For a meaningful alignment we need to select corresponding atoms, 
        # which is complex without sequence alignment. 
        # For this step, we will just acknowledge we loaded both.
        print("Ready for structural comparison.")

    except Exception as e:
        print(f"Error handling 1IRK: {e}")

if __name__ == "__main__":
    analyze_structure()
