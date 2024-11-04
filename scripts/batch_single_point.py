import subprocess
from pathlib import Path

jobname = "id_123"

# Create directory structure
go_dir = Path(f"data/{jobname}_go")
base_dir = Path(f"data/{jobname}_sp")
base_dir.mkdir(exist_ok=True)

# molecules is a list of all dir names in base_dir
molecules = [
    {"name": sub_dir.name} for sub_dir in go_dir.iterdir() if sub_dir.is_dir()
]

# Generate calculation commands
for mol in molecules:
    mol_dir = base_dir / mol["name"]
    mol_dir.mkdir(exist_ok=True)

    if mol['name'].endswith("_cation"):
        charge = 1
        multiplicity = 1
    else:
        charge = 0
        multiplicity = 1

    # Single point with higher level method
    sp_cmd = [
        "ophysia", "create",
        "--xyz-file", str(go_dir / mol["name"] / f"molecule.xyz"),
        "--output-dir", str(mol_dir),
        "--functional", "PBE0",  # More accurate functional with dispersion
        "--basis", "def2-TZVPP",    # Larger basis set for final energy
        "--dispersion", "D4",       # Dispersion
        "--calc-type", "SP-CUSTOM",
        "--charge", str(charge),
        "--multiplicity", str(multiplicity)
    ]
    
    print(f"\nGenerating inputs for {mol['name']}:")
    print(" ".join(sp_cmd))
    subprocess.run(sp_cmd)