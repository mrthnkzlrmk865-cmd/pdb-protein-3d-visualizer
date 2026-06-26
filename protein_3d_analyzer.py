"""
Computational Biology: PDB Data Retrieval & 3D Structural Analysis
Author: Murathan Kizilirmak
"""

import os
import urllib.request
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser

def fetch_pdb_structure(pdb_id):
    """
    Dynamically fetches the verified 3D structural PDB file 
    from the official RCSB Protein Data Bank.
    """
    # Standard URL pattern for RCSB PDB files
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    filename = f"{pdb_id.upper()}.pdb"
    
    print(f"[API Request] Connecting to RCSB Protein Data Bank for PDB ID: {pdb_id}...")
    
    try:
        if not os.path.exists(filename):
            urllib.request.urlretrieve(url, filename)
            print(f"-> Successfully downloaded 3D structure: {filename}")
        else:
            print(f"-> Local cache found for {filename}. Skipping download.")
        return filename
    except Exception as e:
        print(f"-> Error fetching data from PDB: {e}")
        return None

def analyze_and_plot_3d_structure(pdb_filename):
    """
    Parses the molecular PDB file using Bio.PDB 
    and renders a 3D stereoscopic structural trajectory of the alpha-carbon backbone.
    """
    if not pdb_filename:
        return
        
    print(f"\n[Bio.PDB Parsing] Extracting atomic coordinates from {pdb_filename}...")
    
    # Initialize Biopython PDB Parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_filename)
    
    # Extract Alpha Carbon (CA) coordinates to map the 3D protein fold backbone
    x_coords = []
    y_coords = []
    z_coords = []
    
    # Iterate through models, chains, residues, and atoms using Biopython structural hierarchy
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:  # CA represents the central Alpha Carbon of an amino acid
                    atom = residue["CA"]
                    coord = atom.get_coord()
                    x_coords.append(coord[0])
                    y_coords.append(coord[1])
                    z_coords.append(coord[2])
                    
    total_residues = len(x_coords)
    
    if total_residues == 0:
        print("-> Warning: No Alpha Carbon coordinates found in this file format.")
        return
        
    print(f"-> Successfully mapped {total_residues} amino acid alpha-carbons in 3D space.")
    
    # Generate 3D Spatial Visualization
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the peptide backbone line (Protein Folding Trajectory)
    ax.plot(x_coords, y_coords, z_coords, color='#1f77b4', linewidth=2.5, label="Peptide Backbone (CA)")
    
    # Scatter plot to emphasize distinct amino acid positions, colored by sequence index
    sc = ax.scatter(x_coords, y_coords, z_coords, c=range(total_residues), 
                    cmap='viridis', s=20, edgecolors='black', alpha=0.8)
    
    # Add stereoscopic visualization details
    ax.set_title(f"3D Protein Folding Visualization (RCSB PDB BackBone Trajectory)\nSource: {pdb_filename}", 
                 fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel("X Spatial Coordinate (Å)", fontsize=10)
    ax.set_ylabel("Y Spatial Coordinate (Å)", fontsize=10)
    ax.set_zlabel("Z Spatial Coordinate (Å)", fontsize=10)
    
    # Colorbar showing sequence flow from N-Terminus to C-Terminus
    cbar = plt.colorbar(sc, ax=ax, pad=0.1)
    cbar.set_label("Amino Acid Sequence Index (N-Terminus -> C-Terminus)", fontsize=10)
    
    ax.legend()
    plt.tight_layout()
    plt.show()

def main():
    print("=" * 70)
    print("   COMPUTATIONAL BIOLOGY: PDB COORD EXTRACTOR & 3D ANALYZER   ")
    print("=" * 70)
    
    # Using PDB ID '1OLG' which is the structural crystal model of human tumor suppressor P53 tetramerization domain
    target_pdb_id = "1olg" 
    
    # Pipeline execution: Download -> Parse -> Analyze -> 3D Render
    pdb_file = fetch_pdb_structure(target_pdb_id)
    analyze_and_plot_3d_structure(pdb_file)
    
    print("\nPipeline execution complete.")
    print("=" * 70)

if __name__ == "__main__":
    main()
