import os
import pandas as pd
from ophysia.output import parse_single_point_energy
from decimal import Decimal
import re

def is_cation(molecule_name):
    """Check if a molecule name contains 'cation'"""
    return 'cation' in molecule_name.lower()

def parse_energies(filepath):
    """
    Extract single point energy, Gibbs free energy and enthalpy from output file.
    Returns a tuple of (single_point_energy, gibbs_energy, enthalpy)
    """
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
            # Pattern for Gibbs free energy
            gibbs_pattern = r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)\s+Eh"
            gibbs_matches = re.findall(gibbs_pattern, content)
            
            # Pattern for enthalpy
            enthalpy_pattern = r"Total enthalpy\s+\.\.\.\s+(-?\d+\.\d+)\s+Eh"
            enthalpy_matches = re.findall(enthalpy_pattern, content)
            
            # Get single point energy using existing function
            sp_energy = parse_single_point_energy(filepath)
            
            gibbs = Decimal(gibbs_matches[-1]) if gibbs_matches else None
            enthalpy = Decimal(enthalpy_matches[-1]) if enthalpy_matches else None
            
            return sp_energy, gibbs, enthalpy
            
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None, None, None

def collect_energies(directory):
    """
    Walk through subdirectories and collect all energy values.
    Returns a dictionary of {directory_name: {"sp": sp_energy, "gibbs": gibbs_energy, "enthalpy": enthalpy}}
    """
    energies = {}
    
    for dirpath, dirnames, filenames in os.walk(directory):
        if dirpath == directory:
            continue
            
        if 'molecule.out' in filenames:
            dir_name = os.path.basename(dirpath)
            filepath = os.path.join(dirpath, 'molecule.out')
            sp_energy, gibbs, enthalpy = parse_energies(filepath)
            
            if sp_energy is not None:
                energies[dir_name] = {
                    "sp": float(sp_energy),
                    "gibbs": float(gibbs) if gibbs else None,
                    "enthalpy": float(enthalpy) if enthalpy else None
                }
                print(f"Found energies for {dir_name}: SP={sp_energy}, G={gibbs}, H={enthalpy}")

    return energies

def process_directory(directory, base_molecule, base_molecule_protonated):
    """Process a single directory and return results DataFrame"""
    print(f"\nProcessing directory: {directory}")
    
    energies = collect_energies(directory)
    
    if base_molecule not in energies:
        print(f"Warning: Neutral base molecule '{base_molecule}' not found in {directory}")
        return None
        
    if base_molecule_protonated not in energies:
        print(f"Warning: Protonated base molecule '{base_molecule_protonated}' not found in {directory}")
        return None
    
    base_energies = energies[base_molecule]
    base_energies_prot = energies[base_molecule_protonated]
    
    pairs = []
    
    # Convert to DataFrame for easier indexing
    df = pd.DataFrame([
        {
            'molecule': mol,
            'single_point_energy': data['sp'],
            'gibbs_energy': data['gibbs'],
            'enthalpy': data['enthalpy']
        }
        for mol, data in energies.items()
    ])
    
    for i in range(len(df)):
        for j in range(i+1, len(df)):
            mol1 = df['molecule'][i]
            mol2 = df['molecule'][j]

            if mol1 == base_molecule or mol2 == base_molecule or \
               mol1 == base_molecule_protonated or mol2 == base_molecule_protonated:
                continue
            
            if is_cation(mol1) != is_cation(mol2):
                # Calculate SP energy differences
                sp_diff = (float(df['single_point_energy'][i]) + 
                         float(df['single_point_energy'][j]) - 
                         base_energies['sp'])
                sp_diff_prot = (float(df['single_point_energy'][i]) + 
                              float(df['single_point_energy'][j]) - 
                              base_energies_prot['sp'])
                
                # Calculate Gibbs energy differences if available
                g_diff = None
                g_diff_prot = None
                if all(x is not None for x in [df['gibbs_energy'][i], df['gibbs_energy'][j], base_energies['gibbs']]):
                    g_diff = (float(df['gibbs_energy'][i]) + 
                            float(df['gibbs_energy'][j]) - 
                            base_energies['gibbs'])
                    g_diff_prot = (float(df['gibbs_energy'][i]) + 
                                 float(df['gibbs_energy'][j]) - 
                                 base_energies_prot['gibbs'])
                
                # Calculate enthalpy differences if available
                h_diff = None
                h_diff_prot = None
                if all(x is not None for x in [df['enthalpy'][i], df['enthalpy'][j], base_energies['enthalpy']]):
                    h_diff = (float(df['enthalpy'][i]) + 
                            float(df['enthalpy'][j]) - 
                            base_energies['enthalpy'])
                    h_diff_prot = (float(df['enthalpy'][i]) + 
                                 float(df['enthalpy'][j]) - 
                                 base_energies_prot['enthalpy'])
                
                # Convert all to kcal/mol
                conv = 627.509
                pairs.append({
                    'pair': f"{mol1} + {mol2}",
                    'energy_diff': sp_diff * conv,
                    'energy_diff_prot': sp_diff_prot * conv,
                    'gibbs_diff': g_diff * conv if g_diff is not None else None,
                    'gibbs_diff_prot': g_diff_prot * conv if g_diff_prot is not None else None,
                    'enthalpy_diff': h_diff * conv if h_diff is not None else None,
                    'enthalpy_diff_prot': h_diff_prot * conv if h_diff_prot is not None else None
                })
    
    return pd.DataFrame(pairs) if pairs else None

def main(base_dir, molecule_prefix, molecule_id, methods_dir_go, methods_dir_sp):
    """Main function with updated results handling"""
    sp_dir = os.path.join(base_dir, methods_dir_sp, f"{molecule_prefix}_sp_only")
    gosp_dir = os.path.join(base_dir, methods_dir_go, f"{molecule_prefix}_go_sp")
    
    print(f"Analyzing directories:")
    print(f"SP directory: {sp_dir}")
    print(f"GO+SP directory: {gosp_dir}")
    
    base_molecule = molecule_id
    base_molecule_protonated = f"{molecule_id}_protonated"
    
    results = {}
    
    if os.path.exists(sp_dir):
        sp_results = process_directory(sp_dir, base_molecule, base_molecule_protonated)
        if sp_results is not None:
            results['sp'] = sp_results
    
    if os.path.exists(gosp_dir):
        gosp_results = process_directory(gosp_dir, base_molecule, base_molecule_protonated)
        if gosp_results is not None:
            results['gosp'] = gosp_results
    
    if results:
        final_df = pd.DataFrame()
        all_pairs = set()
        for df in results.values():
            all_pairs.update(df['pair'])
        
        final_df['pair'] = list(all_pairs)
        
        # Add SP results
        if 'sp' in results:
            sp_df = results['sp']
            final_df['BDE_SP*'] = final_df['pair'].map(
                sp_df.set_index('pair')['energy_diff']).round(1)
            final_df['BDE_SP*(prot)'] = final_df['pair'].map(
                sp_df.set_index('pair')['energy_diff_prot']).round(1)
            # Add Gibbs and enthalpy for SP
            # final_df['G_SP*'] = final_df['pair'].map(
            #     sp_df.set_index('pair')['gibbs_diff']).round(1)
            # final_df['G_SP*(prot)'] = final_df['pair'].map(
            #     sp_df.set_index('pair')['gibbs_diff_prot']).round(1)
            # final_df['H_SP*'] = final_df['pair'].map(
            #     sp_df.set_index('pair')['enthalpy_diff']).round(1)
            # final_df['H_SP*(prot)'] = final_df['pair'].map(
            #     sp_df.set_index('pair')['enthalpy_diff_prot']).round(1)
        
        # Add GO+SP results
        if 'gosp' in results:
            gosp_df = results['gosp']
            final_df['BDE_GO&SP*'] = final_df['pair'].map(
                gosp_df.set_index('pair')['energy_diff']).round(1)
            final_df['BDE_GO&SP*(prot)'] = final_df['pair'].map(
                gosp_df.set_index('pair')['energy_diff_prot']).round(1)
            # Add Gibbs and enthalpy for GO+SP
            final_df['G_GO&SP*'] = final_df['pair'].map(
                gosp_df.set_index('pair')['gibbs_diff']).round(1)
            final_df['G_GO&SP*(prot)'] = final_df['pair'].map(
                gosp_df.set_index('pair')['gibbs_diff_prot']).round(1)
            final_df['H_GO&SP*'] = final_df['pair'].map(
                gosp_df.set_index('pair')['enthalpy_diff']).round(1)
            final_df['H_GO&SP*(prot)'] = final_df['pair'].map(
                gosp_df.set_index('pair')['enthalpy_diff_prot']).round(1)
        
        output_file = os.path.join(base_dir, f"{molecule_prefix}_energy_results.csv")
        final_df.to_csv(output_file, index=False)
        
        print("\nFinal Results:")
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