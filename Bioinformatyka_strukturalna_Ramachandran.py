import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB
from Bio.PDB.Polypeptide import three_to_one
import sys

# funkcja obliczajaca katy phi i psi dla kazdego aminokwasu w strukturze
def calculate_phi_psi(structure):
    angles = {'phi': [], 'psi': []}
    
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)# buduje polipeptydy z lancucha w strukturze bia³ka
            
            for poly_index, poly in enumerate(polypeptides):
                phi_psi_angles = poly.get_phi_psi_list()# wez phi i psi z residue
                
                for res_index, angles_tuple in enumerate(phi_psi_angles):
                    phi, psi = angles_tuple
                    if phi is not None and psi is not None:# sprawdzanie czy oba phi i psi istnieja w residue
                        angles['phi'].append(np.degrees(phi))
                        angles['psi'].append(np.degrees(psi))
                        
    return angles

# funkcja bioraca katy i strukture drugorzedowa i robiaca z tego ramachadran plot
def plot_ramachandran(angles, sec_structure, save_path=None):
    colors = {'H': 'red', 'E': 'blue', 'C': 'green'}  # Helix, Sheet, Coil
                                                      # Inicjalizacja slownika plotu
    
    plt.figure(figsize=(10, 8))# inicjalizacja figury
    
    for sec, color in colors.items():# przechodzi przez wszystkie struktury ze slownika colors, ekstraktuje indeksy dla kazdej z nich i tworzy scatterplota 
        indices = [i for i, sec_type in enumerate(sec_structure) if sec_type == sec]
        plt.scatter(np.array(angles['phi'])[indices], np.array(angles['psi'])[indices], label=sec, color=color)

    # dodaje cechy plota
    plt.title('Ramachandran Plot')
    plt.xlabel('Phi (degrees)')
    plt.ylabel('Psi (degrees)')
    plt.legend()
    
    # zapisz lub pokaz plot
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
    
    secondary_structure = ['C'] * len(phi_psi_angles['phi'])
    
    plot_ramachandran(phi_psi_angles, secondary_structure)
