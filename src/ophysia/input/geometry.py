from pathlib import Path
from typing import List, Optional

class GeometryBuilder:
    def __init__(self, xyz_file: Path):
        self.xyz_file = xyz_file
        
    def create_input_geometry(self, charge: int = 0, multiplicity: int = 1) -> List[str]:
        """Creates ORCA-formatted geometry section from XYZ file.

        Args:
            charge (int): Charge of molecule
            multiplicity (int): Multiplicity of molecule
        
        Returns:
            List[str]: List of ORCA-formatted geometry lines
        """
        
        geometry_lines = []
        geometry_lines.append(f"* xyz {charge} {multiplicity}")
        
        with open(self.xyz_file) as f:
            # Skip first two lines of XYZ file
            next(f)
            next(f)
            for line in f:
                geometry_lines.append(line.strip())
                
        geometry_lines.append("*")
        return geometry_lines