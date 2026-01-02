
import os
import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser

def find_active_site_cluster(cif_path):
    parser = MMCIFParser()
    structure = parser.get_structure("Toxin", cif_path)
    
    residues = [r for r in structure.get_residues()]
    
    his_residues = [r for r in residues if r.get_resname() == "HIS"]
    acidic_residues = [r for r in residues if r.get_resname() in ["ASP", "GLU"]]
    
    print(f"Scanning for catalytic triad/dyad candidates in {os.path.basename(cif_path)}...")
    print(f"Total HIS: {len(his_residues)}")
    print(f"Total ASP/GLU: {len(acidic_residues)}")
    
    possible_sites = []
    
    # Simple distance check: Find His within 4-5A of Asp/Glu
    for h in his_residues:
        h_atoms = [a for a in h if a.element != "H"] # Skip hydrogens
        if not h_atoms: continue
        h_center = np.mean([a.get_vector().get_array() for a in h_atoms], axis=0)
        
        partners = []
        for a in acidic_residues:
            a_atoms = [at for at in a if at.element != "H"]
            if not a_atoms: continue
            a_center = np.mean([at.get_vector().get_array() for at in a_atoms], axis=0)
            
            dist = np.linalg.norm(h_center - a_center)
            if dist < 6.0: # Generous threshold for side chain interaction
                partners.append((a, dist))
        
        if partners:
            possible_sites.append((h, partners))

    # Sort by number of partners (triads often have His interacting with Asp/Glu)
    possible_sites.sort(key=lambda x: len(x[1]), reverse=True)
    
    print("\n--- Potential Active Site Clusters ---")
    for his, partner_list in possible_sites[:5]: # Show top 5
        print(f"Histidine {his.get_resname()}{his.get_id()[1]} interacts with:")
        for p, d in partner_list:
            print(f"  - {p.get_resname()}{p.get_id()[1]} (Dist: {d:.2f} A)")
            
    return possible_sites

if __name__ == "__main__":
    cif_file = r"c:\Users\nanos\OneDrive\Documentos\Curitas\fold_rinconsito_1_model_0.cif"
    find_active_site_cluster(cif_file)
