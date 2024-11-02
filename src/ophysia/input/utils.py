import subprocess
from pathlib import Path
from typing import Optional
from rdkit import Chem

def is_obabel_available() -> bool:
    """Check if OpenBabel is available on the system."""
    try:
        subprocess.run(["obabel", "--help"], 
                      stdout=subprocess.DEVNULL, 
                      stderr=subprocess.DEVNULL)
        return True
    except FileNotFoundError:
        return False

def validate_smiles(smiles: str) -> Optional[str]:
    """Validate SMILES string using RDKit.
    
    Args:
        smiles (str): SMILES string to validate
    
    Returns:
        Optional[str]: Canonical SMILES if valid, None if invalid
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    return Chem.MolToSmiles(mol, canonical=True)

def smiles_to_xyz(smiles: str, output_file: Path) -> Path:
    """Convert SMILES to 3D XYZ using OpenBabel.

    Args:
        smiles (str): SMILES string of molecule
        output_file (Path): Path to write XYZ file
    
    Returns:
        output_file (Path): Path to written XYZ file

    Raises:
        ValueError: If SMILES string is invalid
        RuntimeError: If conversion fails
    """
    canonical_smiles = validate_smiles(smiles)
    if canonical_smiles is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    try:
        # Run OpenBabel with output capture
        result = subprocess.run(
            [
                "obabel",
                "-:"+canonical_smiles,  # Use canonical SMILES
                "-oxyz",
                "--gen3d",
                "-O", str(output_file)
            ],
            check=True,
            capture_output=True,
            text=True
        )
        
        # Check if output file was created and is not empty
        if not output_file.exists() or output_file.stat().st_size == 0:
            raise RuntimeError("OpenBabel failed to generate XYZ file")
            
        # Verify the output file contains valid XYZ format
        content = output_file.read_text().strip().split('\n')
        if len(content) < 3:  # XYZ files must have at least 3 lines
            raise RuntimeError("Generated XYZ file is invalid")
            
        try:
            natoms = int(content[0])  # First line should be number of atoms
            if len(content) != natoms + 2:  # Check if number of lines matches
                raise RuntimeError("Generated XYZ file has incorrect format")
        except ValueError:
            raise RuntimeError("Generated XYZ file has invalid atom count")
            
        return output_file
        
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to convert SMILES to XYZ: {e.stderr}")
    except Exception as e:
        raise RuntimeError(f"Unexpected error during SMILES conversion: {str(e)}")
