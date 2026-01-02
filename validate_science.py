
import numpy as np
import os
from Bio.PDB import PDBParser, NeighborSearch, PDBIO, SASA, Structure, Model, Chain
from Bio import PDB

def validate_science():
    print("--- Scientific Validation Protocol ---")
    
    receptor_file = "pdb1irk.ent"
    ligand_file = "toxin_aligned_to_1irk.pdb"
    
    if not os.path.exists(ligand_file):
        print("Error: Ligand file missing.")
        return

    parser = PDBParser(QUIET=True)
    receptor_struct = parser.get_structure("Receptor", receptor_file)
    ligand_struct = parser.get_structure("Ligand", ligand_file)
    
    # 1. Reproduce Docking Pose
    print("[1] Reproducing Optimized Docking Pose...")
    rec_atoms = [a for a in receptor_struct.get_atoms() if a.name == 'CA']
    lig_atoms = [a for a in ligand_struct.get_atoms() if a.name == 'CA']
    
    ns = NeighborSearch(rec_atoms)
    
    # Simple De-clash & Refine logic (same as before to ensure consistency)
    current_coords = [a.coord.copy() for a in lig_atoms]
    rec_center = np.mean([a.coord for a in rec_atoms], axis=0)
    lig_center = np.mean([a.coord for a in lig_atoms], axis=0)
    col_vec = (lig_center - rec_center) / np.linalg.norm(lig_center - rec_center)

    # De-clash
    for _ in range(20):
        clashes = 0
        for coord in current_coords:
            if ns.search(coord, 4.0, level='A'): collisions = True; break
        # Check count
        count = 0
        for c in current_coords:
            if [n for n in ns.search(c, 4.0, level='A')]: count+=1
        if count == 0: break
        current_coords = [c + (col_vec * 1.5) for c in current_coords]
        
    # Refine
    best_contacts = 0
    for _ in range(50):
        move = (np.random.rand(3) - 0.5) * 1.0
        trial = [c + move for c in current_coords]
        # Score
        contacts = 0
        clash = False
        for c in trial:
            in_contact = len([n for n in ns.search(c, 7.0, level='A')]) > 0
            in_clash = len([n for n in ns.search(c, 4.0, level='A')]) > 0
            if in_clash: clash = True; break
            if in_contact: contacts += 1
        
        if not clash: 
            # Simple greedy maximization of contacts
            if contacts > best_contacts or (contacts == best_contacts and np.random.rand() > 0.5):
                best_contacts = contacts
                current_coords = trial

    print("    Docking Optimization Complete.")
    
    # Update Ligand Structure in memory with new coordinates
    # Note: We are updating atoms in ligand_struct directly
    # Need to update ALL atoms, but we only simulated CA. 
    # For BSA and Contacts, we ideally need all atoms. 
    # Valid approximation: Apply the translation vector of the Center of Mass to the whole structure.
    
    original_ca_center = np.mean([a.coord for a in lig_atoms if a.name=='CA'], axis=0)
    new_ca_center = np.mean(current_coords, axis=0)
    translation = new_ca_center - original_ca_center
    
    print(f"    Applying Translation: {translation}")
    
    for atom in ligand_struct.get_atoms():
        atom.set_coord(atom.coord + translation)
        
    # Save Final Docked
    io = PDBIO()
    io.set_structure(ligand_struct)
    io.save("toxin_final_docked_validated.pdb")
    
    # 2. Identify Interactions (<= 4.0 A)
    # Re-build NS with all receptor atoms for accuracy
    print("\n[2] Analyzing Real Interface Contacts (< 4.0 A)...")
    rec_all_atoms = list(receptor_struct.get_atoms())
    ns_full = NeighborSearch(rec_all_atoms)
    
    contacts_found = []
    
    toxin_residues = list(ligand_struct.get_residues())
    
    for res in toxin_residues:
        for atom in res:
            neighbors = ns_full.search(atom.coord, 4.0, level='R') # Search for Residues
            for n_res in neighbors:
                # Avoid duplicates
                pair = (f"{res.get_resname()}{res.get_id()[1]}", f"{n_res.get_resname()}{n_res.get_id()[1]}")
                if pair not in contacts_found:
                    contacts_found.append(pair)
    
    print(f"    Total Contact Pairs: {len(contacts_found)}")
    print("    Sample Contacts:")
    for t, r in contacts_found[:10]:
        print(f"      Toxin {t} -- 1IRK {r}")
        
    # 3. BSA Calculation (Using SASA)
    print("\n[3] Calculating Buried Surface Area (BSA)...")
    sr = SASA.ShrakeRupley()
    
    # SASA Complex
    # Create a combined structure object (hacky but functional for BSA)
    # We can just sum atoms? No, SASA needs structure hierarchy.
    # We will compute SASA separately: 
    # BSA = (SASA_A + SASA_B - SASA_Complex) / 2
    # But SASA_Complex requires them in same structure object to detect occlusion.
    
    # Simply measuring SASA of Ligand in isolation vs Ligand in presence of Receptor neighbors?
    # Let's try to calculate SASA for ligand, then SASA for ligand with receptor atoms added to calculation?
    # Biopython SASA works on entities.
    
    # Estimate BSA simply: Count interface atoms * Avg Area subtraction?
    # More robust: 
    # Create dummy structure with both chains.
    
    class DummyStructure:
        def __init__(self, atoms):
            self.atoms = atoms
        def get_atoms(self):
            return self.atoms
        def __iter__(self):
            return iter(self.atoms)

    # Combined atoms list
    all_atoms = list(ligand_struct.get_atoms()) + list(receptor_struct.get_atoms())
    
    # SR computes on structure. Let's try to run it on the separate structures first to establish baseline.
    try:
        sr.compute(ligand_struct, level='S') 
        sasa_lig = ligand_struct.sasa
        
        sr.compute(receptor_struct, level='S')
        sasa_rec = receptor_struct.sasa

        # For complex, we need to trick it or build a proper structure. 
        structure_complex = Structure.Structure("Complex")
        model_c = Model.Model(0)
        structure_complex.add(model_c)
        # Copy chains
        # Deep copy needed to avoid modifying originals? 
        # Python copy module
        import copy
        chain_l = copy.deepcopy(list(ligand_struct.get_chains())[0])
        chain_l.id = 'L'
        model_c.add(chain_l)
        
        # Receptor often has many chains. Add first relevant one (Kinase domain usually A or B)
        # 1IRK has chain A as kinase.
        chain_r = copy.deepcopy(list(receptor_struct.get_chains())[0])
        chain_r.id = 'R'
        model_c.add(chain_r)
        
        sr.compute(structure_complex, level='S')
        sasa_complex = structure_complex.sasa
        
        bsa = (sasa_lig + sasa_rec - sasa_complex) / 2
        print(f"    SASA Ligand: {sasa_lig:.2f} A^2")
        print(f"    SASA Receptor: {sasa_rec:.2f} A^2")
        print(f"    SASA Complex: {sasa_complex:.2f} A^2")
        print(f"    Calculated BSA: {bsa:.2f} A^2")
        
    except Exception as e:
        print(f"    BSA Calculation Error: {e}")
        bsa = 0.0

    # 4. Mutation Check (Clashes)
    print("\n[4] Checking Mutations H38A / E58A...")
    # Check steric clashes for hypothetical Ala at 38 and 58
    # Ala is smaller, so if the Wild Type doesn't clash severely, Ala won't.
    # We check the WT residues 38 and 58 for clashes < 2.5 A (severe).
    
    targets = [38, 58]
    clashes_found = False
    for res in toxin_residues:
        if res.get_id()[1] in targets:
            print(f"    Checking {res.get_resname()}{res.get_id()[1]}...")
            for atom in res:
                # Check distance to Receptor
                neighbors = ns_full.search(atom.coord, 2.5, level='A') # Strict clash
                if neighbors:
                    print(f"      CLASH ALERT: {atom.name} hits Receptor!")
                    clashes_found = True
    
    if not clashes_found:
        print("    No severe steric clashes found for WT residues.")
        print("    Implies H38A/E58A (smaller side chains) are structurally safe.")
    else:
        print("    Warning: WT residues clash. Mutations might actually RELIEVE this.")

if __name__ == "__main__":
    validate_science()
