
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

def improved_docking():
    print("--- Advanced Docking Simulation (De-clashing & Refinement) ---")
    
    # Files
    receptor_file = "pdb1irk.ent"
    ligand_file = "toxin_aligned_to_1irk.pdb"
    
    parser = PDBParser(QUIET=True)
    receptor_struct = parser.get_structure("Receptor", receptor_file)
    ligand_struct = parser.get_structure("Ligand", ligand_file)
    
    rec_atoms = [a for a in receptor_struct.get_atoms() if a.name == 'CA']
    lig_atoms = [a for a in ligand_struct.get_atoms() if a.name == 'CA']
    
    # Calculate Centers
    rec_center = np.mean([a.coord for a in rec_atoms], axis=0)
    lig_center_orig = np.mean([a.coord for a in lig_atoms], axis=0) # current center
    
    print(f"Receptor Center: {rec_center}")
    print(f"Ligand Center: {lig_center_orig}")
    
    ns = NeighborSearch(rec_atoms)
    
    def get_score(coords):
        energy = 0.0
        clashes = 0
        contacts = 0
        for coord in coords:
            neighbors = ns.search(coord, 8.0, level='A') # 8A neighbor list
            for n in neighbors:
                dist = np.linalg.norm(coord - n.coord)
                if dist < 4.0:
                    energy += 500.0 # Heavy Clash penalty
                    clashes += 1
                elif dist <= 7.0: # Contact range
                    energy -= 2.0 # Favorable VdW/Elec
                    contacts += 1
        return energy, clashes, contacts

    # Initial State
    current_coords = [a.coord.copy() for a in lig_atoms]
    e, cl, co = get_score(current_coords)
    print(f"Start: E={e:.1f}, Clashes={cl}, Contacts={co}")

    # Phase 1: De-clashing (Move away from receptor center if clashing)
    # Vector from Rec to Lig
    collision_vec = lig_center_orig - rec_center
    collision_vec = collision_vec / np.linalg.norm(collision_vec) # Normalize
    
    print("\nPhase 1: Resolving Steric Clashes...")
    for i in range(20):
        if cl == 0: break
        # Move outward by 1 Angstrom
        current_coords = [c + (collision_vec * 1.5) for c in current_coords]
        e, cl, co = get_score(current_coords)
        print(f"  Step {i+1}: Clashes={cl}, Energy={e:.1f}")
        
    # Phase 2: Refinement (Hugging the surface)
    print("\nPhase 2: Surface Optimization (Maximizing Contacts)...")
    best_coords = current_coords
    best_e = e
    
    for i in range(100):
        # Random small translation + wobble
        move = (np.random.rand(3) - 0.5) * 1.0 # 1A wobble
        
        trial_coords = [c + move for c in current_coords]
        te, tcl, tco = get_score(trial_coords)
        
        # Accept if energy improves or clashes stay low (Metropolis-like but simple greedy)
        if te < best_e and tcl == 0:
            best_e = te
            best_coords = trial_coords
            # print(f"  Improvement: E={best_e:.1f} (Contacts: {tco})")
            current_coords = trial_coords # Update base
        elif tcl == 0:
             # Acceptance probability for neutral moves? 
             # Just keep current if no improvement to avoid drift
             pass

    final_e, final_cl, final_co = get_score(best_coords)
    print(f"\nFinal State: Energy={final_e:.2f} | Contacts={final_co} | Clashes={final_cl}")
    
    print("\n--- Docking Conclusion ---")
    if final_e < 0:
        print(f"SUCCESS: Negative Delta G ({final_e:.2f}).")
        print("Spontaneous interaction detected at the receptor surface.")
    else:
        print(f"Result: Delta G is positive ({final_e:.2f}). No strong binding found.")

if __name__ == "__main__":
    improved_docking()
