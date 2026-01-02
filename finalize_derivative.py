
import numpy as np
import os
from Bio import PDB
from Bio.PDB import PDBParser, NeighborSearch, PDBIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1

def finalize_derivative():
    print("--- Finalizing Loxosceles-Insulin-1 Derivative ---")
    
    # 1. Re-Run Docking to get Final Coordinates
    receptor_file = "pdb1irk.ent"
    ligand_file = "toxin_aligned_to_1irk.pdb"
    
    if not os.path.exists(ligand_file): 
        print("Error: Aligned toxin file missing. Run alignment first.")
        return

    parser = PDBParser(QUIET=True)
    receptor_struct = parser.get_structure("Receptor", receptor_file)
    ligand_struct = parser.get_structure("Ligand", ligand_file)
    
    rec_atoms = [a for a in receptor_struct.get_atoms() if a.name == 'CA']
    lig_atoms = [a for a in ligand_struct.get_atoms() if a.name == 'CA']
    
    # Helper for simple scoring/docking (Same as before)
    ns = NeighborSearch(rec_atoms)
    
    def get_score_and_contacts(coords):
        energy = 0.0
        clashes = 0
        contacts = set() # Store indices of contacting ligand atoms
        
        for idx, coord in enumerate(coords):
            neighbors = ns.search(coord, 8.0, level='A') 
            for n in neighbors:
                dist = np.linalg.norm(coord - n.coord)
                if dist < 4.0:
                    energy += 500.0
                    clashes += 1
                elif dist <= 7.0:
                    energy -= 2.0
                    contacts.add(idx) # Ligand atom index
        return energy, clashes, contacts

    # Initial Coords
    current_coords = [a.coord.copy() for a in lig_atoms]
    rec_center = np.mean([a.coord for a in rec_atoms], axis=0)
    lig_center = np.mean([a.coord for a in lig_atoms], axis=0)
    
    # Fast De-clash
    col_vec = (lig_center - rec_center) / np.linalg.norm(lig_center - rec_center)
    for _ in range(20):
        e, cl, _ = get_score_and_contacts(current_coords)
        if cl == 0: break
        current_coords = [c + (col_vec * 1.5) for c in current_coords]
        
    # Fast Surface Opt
    best_e = e
    for _ in range(100):
        move = (np.random.rand(3) - 0.5) * 1.0
        trial_coords = [c + move for c in current_coords]
        te, tcl, _ = get_score_and_contacts(trial_coords)
        if te < best_e and tcl == 0:
            best_e = te
            current_coords = trial_coords
            
    print(f"Docking Complete. Final Delta G: {best_e:.2f}")
    
    # Apply to Structure & Save
    for i, atom in enumerate(lig_atoms):
        atom.set_coord(current_coords[i])
        # Note: This only moves CA atoms. For a full PDB save we'd need to move residues. 
        # But for sequence analysis, we just need the IDs which we have.
        # For simplicity in this script, we trust the IDs and the CA positions for interface check.
    
    # 2. Analyze Interface
    _, _, contact_indices = get_score_and_contacts(current_coords)
    contact_residues = []
    
    print("\n[Interface Residues - Toxin]")
    # Get residue objects corresponding to atoms
    lig_residues = list(ligand_struct.get_residues()) 
    # Note: lig_atoms list was built from get_atoms(), order usually preserved but let's be safe
    # Actually, we filtered for CA. let's assume 1:1 map for 304 atoms -> 304 residues
    
    # Map index back to residue info
    for idx in contact_indices:
        # PDB iteration is deterministic
        res = lig_residues[idx]
        res_id = res.get_id()[1]
        res_name = res.get_resname()
        contact_residues.append((res_name, res_id))
        # print(f"  - {res_name} {res_id}")
        
    print(f"Total interacting residues: {len(contact_residues)}")
    
    # Check Active Site presence
    active_site_ids = [38, 58]
    for rname, rid in contact_residues:
        if rid in active_site_ids:
            print(f"  >>> CRITICAL: Active Site Residue {rname}{rid} is at the interface!")

    # 3. Generate Final Sequence (Mutagenesis)
    print("\n[Mutagenesis Design]")
    print("Target: Deactivate Sphingomyelinase D activity.")
    print("Strategy: Point mutations H38A and E58A.")
    
    # Extract Sequence
    ppb = PDB.PPBuilder()
    seq_records = ppb.build_peptides(ligand_struct)
    if not seq_records:
        print("Error: Could not extract sequence.")
        return
        
    original_seq = seq_records[0].get_sequence() # Fixed attribute access
    mutable_seq = list(original_seq)
    
    # Apply Mutations (0-indexed list vs 1-indexed PDB)
    # BEWARE: PDB residues might not start at 1 or might have gaps.
    # The 'toxin' structure is from AlphaFold/ESMFold likely clean 1-indexed.
    # Let's verify by checking residue index.
    
    def mutate(res_num, target_char):
        # Find index in list that matches PDB ID
        found = False
        for i, res in enumerate(lig_residues):
            if res.get_id()[1] == res_num:
                old = mutable_seq[i]
                mutable_seq[i] = target_char
                print(f"  Mutation Applied: {res.get_resname()}{res_num} -> {target_char} (Was {old})")
                found = True
                break
        if not found:
            print(f"  Warning: Residue {res_num} not found for mutation.")

    mutate(38, 'A') # HIS -> ALA
    mutate(58, 'A') # GLU -> ALA
    
    final_seq_str = "".join(mutable_seq)
    
    print("\n[Derivado Loxosceles-Insulina-1 Sequence]")
    print(">" + "Derivado_Loxosceles_Insulina_1_H38A_E58A")
    print(final_seq_str)
    
    # Write to FASTA
    with open("derivado_final.fasta", "w") as f:
        f.write(">Derivado_Loxosceles_Insulina_1_H38A_E58A\n")
        f.write(final_seq_str + "\n")
    print("\nSaved to 'derivado_final.fasta'.")

if __name__ == "__main__":
    finalize_derivative()
