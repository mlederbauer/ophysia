import os
import pandas as pd
from ophysia.output import parse_single_point_energy

def collect_energies(root_dir):
    """
    Walk through subdirectories and collect energy values.
    Returns a dictionary of {directory_name: energy_value}
    """
    energies = {}
    
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
                energies[dir_name] = energy

    return energies

def main(root_dir):
    # Collect all energies
    energies = collect_energies(root_dir)
    
    # Convert to DataFrame
    df = pd.DataFrame(
        [(mol, energy) for mol, energy in energies.items()],
        columns=['molecule', 'single_point_energy']
    )
    
    # Save to CSV
    output_file = f'{root_dir}/single_point_energies.csv'
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    print("\nFirst few entries:")
    print(df.head())

    # generate new df that is the difference row 0 - row 1, column 0 - column 1, 0 - 2, etc.
    df_new = pd.DataFrame()
    for i in range(len(df)):
        for j in range(i+1, len(df)):
            df_new.loc[f"{df['molecule'][i]} + {df['molecule'][j]}", 'energy_diff'] = df['single_point_energy'][i] + df['single_point_energy'][j]
    
    print(df_new)

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory_path>")
        sys.exit(1)
    
    root_dir = sys.argv[1]
    if not os.path.isdir(root_dir):
        print(f"Error: {root_dir} is not a valid directory")
        sys.exit(1)
        
    main(root_dir)