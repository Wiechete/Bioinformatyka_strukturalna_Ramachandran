import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB
from Bio.PDB.Polypeptide import three_to_one
import sys

def calculate_phi_psi(structure):
    angles = {'phi': [], 'psi': []}
    
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            
            for poly_index, poly in enumerate(polypeptides):
                phi_psi_angles = poly.get_phi_psi_list()
                
                for res_index, angles_tuple in enumerate(phi_psi_angles):
                    phi, psi = angles_tuple
                    if phi is not None and psi is not None:
                        angles['phi'].append(np.degrees(phi))
                        angles['psi'].append(np.degrees(psi))
                        
    return angles

def plot_ramachandran(angles, sec_structure, save_path=None):
    colors = {'H': 'red', 'E': 'blue', 'C': 'green'}  # Helix, Sheet, Coil
    
    plt.figure(figsize=(10, 8))
    
    for sec, color in colors.items():
        indices = [i for i, sec_type in enumerate(sec_structure) if sec_type == sec]
        plt.scatter(np.array(angles['phi'])[indices], np.array(angles['psi'])[indices], label=sec, color=color)

    plt.title('Ramachandran Plot')
    plt.xlabel('Phi (degrees)')
    plt.ylabel('Psi (degrees)')
    plt.legend()
    
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    filename = os.path.basename(pdb_file)
    pdb_id = filename[0:4]
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)
    
    phi_psi_angles = calculate_phi_psi(structure)
    
    # Assuming you have information about secondary structure (replace with your own method)
    # Example: Assign 'H' for helix, 'E' for sheet, 'C' for coil
    secondary_structure = ['C'] * len(phi_psi_angles['phi'])
    
    plot_ramachandran(phi_psi_angles, secondary_structure)
