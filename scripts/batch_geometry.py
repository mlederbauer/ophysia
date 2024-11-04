import subprocess
from pathlib import Path
import json
import sys

def load_molecule_config(base_dir):
    """Load molecule configuration from config.json in the base directory."""
    config_path = base_dir / "config.json"
    try:
        with open(config_path) as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: Configuration file not found at {config_path}")
        print("Please ensure config.json exists with molecule definitions.")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in {config_path}")
        sys.exit(1)

def validate_molecule_config(molecules):
    """Validate that all required fields are present in the molecule configurations."""
    required_fields = {"name", "smiles", "charge", "multiplicity"}
    
    for i, mol in enumerate(molecules):
        missing_fields = required_fields - set(mol.keys())
        if missing_fields:
            print(f"Error: Molecule at index {i} is missing required fields: {missing_fields}")
            sys.exit(1)
        
        # Validate data types
        if not isinstance(mol["name"], str):
            print(f"Error: name must be a string for molecule {mol['name']}")
            sys.exit(1)
        if not isinstance(mol["smiles"], str):
            print(f"Error: smiles must be a string for molecule {mol['name']}")
            sys.exit(1)
        try:
            mol["charge"] = int(mol["charge"])
            mol["multiplicity"] = int(mol["multiplicity"])
        except (ValueError, TypeError):
            print(f"Error: charge and multiplicity must be integers for molecule {mol['name']}")
            sys.exit(1)

def generate_calculations(base_dir, molecules):
    """Generate calculation inputs for each molecule."""
    for mol in molecules:
        mol_dir = base_dir / mol["name"]
        mol_dir.mkdir(exist_ok=True)
        
        # Geometry optimization
        opt_cmd = [
            "ophysia", "create",
            "--smiles", mol["smiles"],
            "--output-dir", str(mol_dir),
            "--functional", "B3LYP",
            "--basis", "def2-SVP",
            "--dispersion", "D4",
            "--calc-type", "OPT-CUSTOM",
            "--charge", str(mol["charge"]),
            "--multiplicity", str(mol["multiplicity"]),
        ]

        print(f"\nGenerating inputs for {mol['name']}:")
        print(" ".join(opt_cmd))
        try:
            subprocess.run(opt_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating inputs for {mol['name']}: {e}")
            continue

def main():
    # Create directory structure
    base_dir = Path("./id_123")
    base_dir.mkdir(exist_ok=True)
    
    # Load and validate configuration
    molecules = load_molecule_config(base_dir)
    validate_molecule_config(molecules)
    
    # Generate calculations
    generate_calculations(base_dir, molecules)

if __name__ == "__main__":
    main()