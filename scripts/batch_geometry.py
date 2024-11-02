import subprocess
from pathlib import Path

# Create directory structure
base_dir = Path("cco_analysis")
base_dir.mkdir(exist_ok=True)

# Define molecules and their properties
molecules = [
    {
        "name": "CCO",
        "smiles": "CCO",
        "charge": 0,
        "multiplicity": 1
    },
    {
        "name": "CC_cation",
        "smiles": "C[CH2]",  # CC fragment
        "charge": 1,
        "multiplicity": 1
    },
    {
        "name": "C_cation",
        "smiles": "[CH3]",  # C fragment
        "charge": 1,
        "multiplicity": 1
    },
    {
        "name": "O",
        "smiles": "O",  # C fragment
        "charge": 0,
        "multiplicity": 1
    },
    {
        "name": "CO",
        "smiles": "CO",  # C fragment
        "charge": 0,
        "multiplicity": 1
    }
]

# Generate calculation commands
for mol in molecules:
    mol_dir = base_dir / mol["name"]
    mol_dir.mkdir(exist_ok=True)
    
    # Geometry optimization
    opt_cmd = [
        "ophysia", "create",
        "--smiles", mol["smiles"],
        "--output-dir", str(mol_dir),
        "--functional", "B3LYP",  # Hybrid functional good for geometry
        "--basis", "def2-SVP",  # Triple-zeta basis set for accuracy
        "--dispersion", "D4",  # D3 dispersion correction
        "--calc-type", "OPT-CUSTOM",
        "--charge", str(mol["charge"]),
        "--multiplicity", str(mol["multiplicity"]),
        # '--n-procs', '4'
    ]

    print(f"\nGenerating inputs for {mol['name']}:")
    print(" ".join(opt_cmd))
    subprocess.run(opt_cmd)
    # print(" ".join(sp_cmd))
    # subprocess.run(sp_cmd)