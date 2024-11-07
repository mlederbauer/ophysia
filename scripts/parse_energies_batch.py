import os
import pandas as pd
from ophysia.output import parse_single_point_energy

def collect_energies(directory):
    """
    Walk through subdirectories and collect energy values.
    Returns a dictionary of {directory_name: energy_value}
    """
    energies = {}
    
    # Walk through all subdirectories
    for dirpath, dirnames, filenames in os.walk(directory):
        # Skip the root directory itself
        if dirpath == directory:
            continue
            
        # Look for molecule.out file
        if 'molecule.out' in filenames:
            # Get the immediate parent directory name
            dir_name = os.path.basename(dirpath)
            
            # Parse the energy from the file
            filepath = os.path.join(dirpath, 'molecule.out')
            energy = parse_single_point_energy(filepath)
            
            if energy is not None:
                energies[dir_name] = float(energy)
                print(f"Found energy for {dir_name}: {energy}")

    return energies

def is_cation(molecule_name):
    """Check if a molecule name contains 'cation'"""
    return 'cation' in molecule_name.lower()

def process_directory(directory, base_molecule, base_molecule_protonated):
    """Process a single directory and return results DataFrame"""
    print(f"\nProcessing directory: {directory}")
    
    # Collect all energies
    energies = collect_energies(directory)
    
    if base_molecule not in energies:
        print(f"Warning: Neutral base molecule '{base_molecule}' not found in {directory}")
        return None
        
    if base_molecule_protonated not in energies:
        print(f"Warning: Protonated base molecule '{base_molecule_protonated}' not found in {directory}")
        return None
    
    base_energy = energies[base_molecule]
    base_energy_protonated = energies[base_molecule_protonated]
    
    # Generate pairs where exactly one molecule contains "cation"
    pairs = []
    
    # Convert to DataFrame for easier indexing
    df = pd.DataFrame(
        [(mol, energy) for mol, energy in energies.items()],
        columns=['molecule', 'single_point_energy']
    )
    
    for i in range(len(df)):
        for j in range(i+1, len(df)):
            mol1 = df['molecule'][i]
            mol2 = df['molecule'][j]

            if mol1 == base_molecule or mol2 == base_molecule or \
               mol1 == base_molecule_protonated or mol2 == base_molecule_protonated:
                continue
            
            # Check if exactly one molecule contains "cation"
            if is_cation(mol1) != is_cation(mol2):  # XOR - one must be cation, one must not be
                # Calculate energy differences relative to both references
                energy_diff = float(df['single_point_energy'][i]) + float(df['single_point_energy'][j]) - base_energy
                energy_diff_prot = float(df['single_point_energy'][i]) + float(df['single_point_energy'][j]) - base_energy_protonated
                
                # Convert to kcal/mol
                energy_diff_kcal = energy_diff * 627.509
                energy_diff_prot_kcal = energy_diff_prot * 627.509
                
                pairs.append({
                    'pair': f"{mol1} + {mol2}",
                    'energy_diff': energy_diff_kcal,
                    'energy_diff_prot': energy_diff_prot_kcal
                })
    
    return pd.DataFrame(pairs) if pairs else None

def main(base_dir, molecule_prefix, molecule_id, methods_dir_go, methods_dir_sp):
    """
    Process both sp_only and go_sp directories for a given molecule
    """
    # Define directories - now directly using the method directories
    sp_dir = os.path.join(base_dir, methods_dir_sp, f"{molecule_prefix}_sp_only")
    gosp_dir = os.path.join(base_dir, methods_dir_go, f"{molecule_prefix}_go_sp")
    
    print(f"Analyzing directories:")
    print(f"SP directory: {sp_dir}")
    print(f"GO+SP directory: {gosp_dir}")
    
    # Define base molecules
    base_molecule = molecule_id
    base_molecule_protonated = f"{molecule_id}_protonated"
    
    results = {}
    
    # Process SP directory if it exists
    if os.path.exists(sp_dir):
        sp_results = process_directory(sp_dir, base_molecule, base_molecule_protonated)
        if sp_results is not None:
            results['sp'] = sp_results
    else:
        print(f"Warning: SP directory not found: {sp_dir}")
    
    # Process GO+SP directory if it exists
    if os.path.exists(gosp_dir):
        gosp_results = process_directory(gosp_dir, base_molecule, base_molecule_protonated)
        if gosp_results is not None:
            results['gosp'] = gosp_results
    else:
        print(f"Warning: GO+SP directory not found: {gosp_dir}")
    
    # Combine results into final DataFrame
    if results:
        final_df = pd.DataFrame()
        
        # Get all unique pairs
        all_pairs = set()
        for df in results.values():
            all_pairs.update(df['pair'])
        
        final_df['pair'] = list(all_pairs)
        
        # Add SP results if available
        if 'sp' in results:
            sp_df = results['sp']
            final_df['BDE_SP*'] = final_df['pair'].map(
                sp_df.set_index('pair')['energy_diff']).round(1)
            final_df['BDE_SP*(prot)'] = final_df['pair'].map(
                sp_df.set_index('pair')['energy_diff_prot']).round(1)
        
        # Add GO+SP results if available
        if 'gosp' in results:
            gosp_df = results['gosp']
            final_df['BDE_GO&SP*'] = final_df['pair'].map(
                gosp_df.set_index('pair')['energy_diff']).round(1)
            final_df['BDE_GO&SP*(prot)'] = final_df['pair'].map(
                gosp_df.set_index('pair')['energy_diff_prot']).round(1)
        
        # Save results
        output_file = os.path.join(base_dir, f"{molecule_prefix}_bde_results.csv")
        final_df.to_csv(output_file, index=False)
        
        print("\nFinal Results:")
        # Sort in descending order by BDE_SP*(prot) column
        final_df = final_df.sort_values(by='BDE_SP*(prot)', ascending=False)
        pd.set_option('display.max_rows', None)
        print(final_df)
        print(f"\nResults saved to: {output_file}")
    else:
        print("\nNo results could be calculated.")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 4:
        print("Usage: python script.py <data_directory> <dir_prefix> <molecule_id>")
        print("Example: python script.py ./data gabaarg gaba")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    dir_prefix = sys.argv[2]    # e.g., 'gabaarg'
    molecule_id = sys.argv[3]   # e.g., 'gaba'
    methods_dir_go = 'GO_PBE0_def2TZVP'
    methods_dir_sp = 'SP_PBE0_def2QZVP'
    
    if not os.path.isdir(base_dir):
        print(f"Error: {base_dir} is not a valid directory")
        sys.exit(1)
        
    main(base_dir, dir_prefix, molecule_id, methods_dir_go, methods_dir_sp)