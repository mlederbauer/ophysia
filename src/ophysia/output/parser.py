import os
import pandas as pd
import re
from decimal import Decimal

def parse_single_point_energy(filepath):
    """Extract single point energy from output file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            # Look for the energy value between "FINAL SINGLE POINT ENERGY" and multiple dashes
            pattern = r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"
            matches = re.findall(pattern, content)
            if matches:
                return Decimal(matches[-1])  # Return the last match
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return None

def parse_energies(filepath):
    """
    Extract Gibbs free energy and enthalpy from output file.
    
    Args:
        filepath (str): Path to the output file
        
    Returns:
        tuple: (free_energy, enthalpy) in Eh units, or (None, None) if parsing fails
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
            
            free_energy = Decimal(gibbs_matches[-1]) if gibbs_matches else None
            enthalpy = Decimal(enthalpy_matches[-1]) if enthalpy_matches else None
            
            return free_energy, enthalpy
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
