import subprocess
from pathlib import Path
import json
import sys
import argparse

def load_molecule_config(config_path):
    """Load molecule configuration from config.json."""
    try:
        with open(config_path) as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: Configuration file not found at {config_path}")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in {config_path}")
        sys.exit(1)

def validate_molecule_config(molecules, input_type):
    """Validate molecule configurations based on input type."""
    if input_type == "smiles":
        required_fields = {"name", "smiles", "multiplicity"}
    else:  # xyz
        required_fields = {"name"}
    
    for i, mol in enumerate(molecules):
        missing_fields = required_fields - set(mol.keys())
        if missing_fields:
            print(f"Error: Molecule at index {i} is missing required fields: {missing_fields}")
            sys.exit(1)
        
        # Validate data types
        if not isinstance(mol["name"], str):
            print(f"Error: name must be a string for molecule {mol['name']}")
            sys.exit(1)
            
        if input_type == "smiles":
            if not isinstance(mol["smiles"], str):
                print(f"Error: smiles must be a string for molecule {mol['name']}")
                sys.exit(1)
            try:
                mol["charge"] = 1 if mol["name"].endswith("_cation") else 0
                mol["multiplicity"] = int(mol["multiplicity"])
            except (ValueError, TypeError):
                print(f"Error: charge and multiplicity must be integers for molecule {mol['name']}")
                sys.exit(1)

def get_calculation_settings(calc_type):
    """Get functional, basis set, and other settings based on calculation type."""
    if calc_type == "GO":
        return {
            "functional": "B3LYP",
            "basis": "def2-SVP",
            "calc_type": "OPT-CUSTOM"
        }
    else:  # SP
        return {
            "functional": "PBE0",
            "basis": "def2-TZVPP",
            "calc_type": "SP-CUSTOM"
        }

def infer_charge_multiplicity(molecule_name):
    """Infer charge and multiplicity from molecule name."""
    if molecule_name.endswith("_cation"):
        return 1, 1
    return 0, 1

def get_molecules_from_directory(xyz_dir):
    """Get molecule list from directory structure."""
    return [{"name": sub_dir.name} for sub_dir in xyz_dir.iterdir() if sub_dir.is_dir()]

def generate_calculations(molecules, input_type, calc_type, base_dir, xyz_dir=None):
    """Generate calculation inputs for each molecule."""
    settings = get_calculation_settings(calc_type)
    
    for mol in molecules:
        mol_dir = base_dir / mol["name"]
        mol_dir.mkdir(exist_ok=True)
        
        # Set charge and multiplicity if not provided
        if "charge" not in mol or "multiplicity" not in mol:
            mol["charge"], mol["multiplicity"] = infer_charge_multiplicity(mol["name"])
        
        # Build common command arguments
        cmd = [
            "ophysia", "create",
            "--output-dir", str(mol_dir),
            "--functional", settings["functional"],
            "--basis", settings["basis"],
            "--dispersion", "D4",
            "--calc-type", settings["calc_type"],
            "--charge", str(mol["charge"]),
            "--multiplicity", str(mol["multiplicity"]),
        ]
        
        # Add input-specific arguments
        if input_type == "smiles":
            cmd.extend(["--smiles", mol["smiles"]])
        else:  # xyz
            xyz_path = xyz_dir / mol["name"] / "molecule.xyz"
            if not xyz_path.exists():
                print(f"Error: XYZ file not found at {xyz_path}")
                continue
            cmd.extend(["--xyz-file", str(xyz_path)])
        
        print(f"\nGenerating inputs for {mol['name']}:")
        print(" ".join(cmd))
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating inputs for {mol['name']}: {e}")
            continue

def main():
    parser = argparse.ArgumentParser(description="Generate quantum chemistry calculations")
    parser.add_argument("--job-id", required=True, help="Job ID (e.g., id_123)")
    parser.add_argument("--calc-type", choices=["GO", "SP"], required=True,
                      help="Calculation type: Geometry Optimization (GO) or Single Point (SP)")
    parser.add_argument("--input-type", choices=["smiles", "xyz"], required=True,
                      help="Input type: SMILES string or XYZ file")
    parser.add_argument("--xyz-dir", help="Path to directory containing XYZ files (required for XYZ input)")
    
    args = parser.parse_args()
    
    # Set up directories
    base_dir = Path(f"data/{args.job_id}")
    base_dir.mkdir(exist_ok=True, parents=True)
    
    # Get molecule configurations
    if args.input_type == "smiles":
        # Check if "data/{args.job_id}/config.json" exists
        if not Path(f"data/{args.job_id}/config.json").exists():
            print(f"Error: Configuration file not found at data/{args.job_id}/config.json")
            raise FileNotFoundError
        molecules = load_molecule_config(f"data/{args.job_id}/config.json")
    else:  # xyz
        if not args.xyz_dir:
            print("Error: XYZ directory is required for XYZ input")
            raise ValueError
        xyz_dir = Path(args.xyz_dir)
        if not xyz_dir.exists():
            print(f"Error: XYZ directory not found at {xyz_dir}")
            raise FileNotFoundError
        molecules = get_molecules_from_directory(xyz_dir)
    
    # Validate configurations
    validate_molecule_config(molecules, args.input_type)
    
    # Generate calculations
    generate_calculations(
        molecules, 
        args.input_type, 
        args.calc_type, 
        base_dir, 
        Path(args.xyz_dir) if args.xyz_dir else None
    )

if __name__ == "__main__":
    main()