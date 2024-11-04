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
            match = re.search(pattern, content)
            if match:
                return Decimal(match.group(1))
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return None
