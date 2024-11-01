import subprocess
from pathlib import Path

def smiles_to_xyz(smiles: str, output_file: Path) -> Path:
    """Convert SMILES to 3D XYZ using OpenBabel.

    Args:
        smiles (str): SMILES string of molecule
        output_file (Path): Path to write XYZ file
    
    Returns:
        output_file (Path): Path to written XYZ file

    Raises:
        RuntimeError: If conversion fails
    """

    try:
        subprocess.run([
            "obabel",
            "-:"+smiles,
            "-oxyz",
            "--gen3d",
            "-O", str(output_file)
        ], check=True)
        return output_file
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to convert SMILES to XYZ: {e}")
