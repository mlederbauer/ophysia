import os
import pandas as pd
from ophysia.output import parse_single_point_energy

def collect_energies(root_dir, base_molecule):
    """
    Walk through subdirectories and collect energy values.
    Returns a dictionary of {directory_name: energy_value}
    """
    energies = {}
    base_energy = None
    
    # Walk through all subdirectories
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Skip the root directory itself
        if dirpath == root_dir:
            continue
            
        # Look for molecule.out file
        if 'molecule.out' in filenames:
            # Get the immediate parent directory name
            dir_name = os.path.basename(dirpath)
            
            # Parse the energy from the file
            filepath = os.path.join(dirpath, 'molecule.out')
            energy = parse_single_point_energy(filepath)
            
            if energy is not None:
                energies[dir_name] = float(energy)  # Convert to float here
                
                # Store base molecule energy separately
                if dir_name == base_molecule:
                    base_energy = float(energy)  # Convert to float here

    return energies, base_energy

def is_cation(molecule_name):
    """Check if a molecule name contains 'cation'"""
    return 'cation' in molecule_name.lower()

def main(root_dir, base_molecule):
    # Collect all energies
    energies, base_energy = collect_energies(root_dir, base_molecule)
    
    if base_energy is None:
        print(f"Warning: Base molecule '{base_molecule}' not found!")
        return
    
    # Convert to DataFrame
    df = pd.DataFrame(
        [(mol, energy) for mol, energy in energies.items()],
        columns=['molecule', 'single_point_energy']
    )
    
    # Save to CSV
    output_file = f'{root_dir}/single_point_energies.csv'
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    print("\nEntries:")
    print(df)

    # Generate new df that contains only pairs where exactly one molecule contains "cation"
    df_new = pd.DataFrame(columns=['pair', 'energy_diff', 'energy_diff_kcalpermol'])
    
    for i in range(len(df)):
        for j in range(i+1, len(df)):
            mol1 = df['molecule'][i]
            mol2 = df['molecule'][j]

            if mol1 == base_molecule or mol2 == base_molecule:
                continue
            
            # Check if exactly one molecule contains "cation"
            if is_cation(mol1) != is_cation(mol2):  # XOR - one must be cation, one must not be
                # Calculate energy difference and subtract base molecule energy
                energy_diff = float(df['single_point_energy'][i]) + float(df['single_point_energy'][j]) - base_energy
                energy_diff_kcalpermol = energy_diff * 627.509  # Convert Hartree to kcal/mol
                df_new.loc[len(df_new)] = [f"{mol1} + {mol2}", energy_diff, energy_diff_kcalpermol]
    
    print("\nPairwise energy differences (with base molecule subtracted):")
    print(df_new)
    
    # Save pairwise results to CSV
    pairwise_output_file = f'{root_dir}/pairwise_energies.csv'
    df_new.to_csv(pairwise_output_file, index=False)
    print(f"\nPairwise results saved to {pairwise_output_file}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python script.py <directory_path> <base_molecule>")
        sys.exit(1)
    
    root_dir = sys.argv[1]
    base_molecule = sys.argv[2]
    
    if not os.path.isdir(root_dir):
        print(f"Error: {root_dir} is not a valid directory")
        sys.exit(1)
        
    main(root_dir, base_molecule)