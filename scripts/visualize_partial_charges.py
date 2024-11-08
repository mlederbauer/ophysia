import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import rdkit
import re
import os

def is_cationic(dirname):
    """Check if molecule is cationic based on directory name"""
    return dirname.endswith(('cation', 'protonated'))

def read_charges(filename):
    """Read Mulliken charges from output file using regex"""
    charges = {}
    try:
        with open(filename, 'r') as f:
            content = f.read()
            # Find the charges section using regex
            pattern = r"MULLIKEN ATOMIC CHARGES\n-+\n([\s\S]+?)\nSum of atomic charges:"
            match = re.search(pattern, content)
            if match:
                charges_text = match.group(1)
                # Parse each charge line
                for line in charges_text.strip().split('\n'):
                    if line.strip():
                        # Match the index and charge value
                        charge_pattern = r"\s*(\d+)\s+[A-Za-z]+\s*:\s*([-+]?\d*\.\d+)"
                        charge_match = re.search(charge_pattern, line)
                        if charge_match:
                            idx = int(charge_match.group(1))
                            charge = float(charge_match.group(2))
                            charges[idx] = charge
    except FileNotFoundError:
        print(f"Error: Could not find file {filename}")
    return charges

def read_xyz(filename):
    """Read XYZ file content"""
    try:
        with open(filename, 'r') as f:
            return f.read()
    except FileNotFoundError:
        print(f"Error: Could not find file {filename}")
        return None

def xyz_to_mol(xyz_content, smiles, is_cation=False):
    """Convert XYZ content to RDKit mol object using SMILES as reference"""
    if not xyz_content:
        return None
    
    # Parse XYZ content
    lines = xyz_content.strip().split('\n')
    n_atoms = int(lines[0])
    atoms = []
    coords = []
    
    for line in lines[2:n_atoms+2]:
        atom, x, y, z = line.split()
        atoms.append(atom)
        coords.append([float(x), float(y), float(z)])
    
    # Create template molecule from SMILES
    template_mol = Chem.MolFromSmiles(smiles)
    if template_mol is None:
        print("Error: Could not parse SMILES")
        return None
    
    # Create mol with explicit bonds from template
    mol = Chem.RWMol()
    
    # Add atoms
    atom_map = {}  # Map from template atom idx to new atom idx
    for i, atom_symbol in enumerate(atoms):
        atom = Chem.rdchem.Atom(atom_symbol)
        new_idx = mol.AddAtom(atom)
        atom_map[i] = new_idx
    
    # Add bonds based on template
    for bond in template_mol.GetBonds():
        begin_idx = atom_map[bond.GetBeginAtomIdx()]
        end_idx = atom_map[bond.GetEndAtomIdx()]
        bond_type = bond.GetBondType()
        mol.AddBond(begin_idx, end_idx, bond_type)
    
    # Set 3D coordinates
    conf = Chem.Conformer(n_atoms)
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf)
    
    # Sanitize and convert to regular mol
    try:
        Chem.SanitizeMol(mol)
        mol = Chem.Mol(mol)
        
        # Generate 2D coordinates for visualization
        AllChem.Compute2DCoords(mol)
        
        return mol
    except:
        print("Error: Could not sanitize molecule. Check atomic coordinates and connectivity.")
        return None

def assign_charges(mol, charges_dict, is_cation=False):
    """Assign partial charges to atoms in the molecule"""
    # Find the most positive H bound to N if cationic
    max_h_charge = -float('inf')
    max_h_idx = None
    
    if is_cation:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'N':
                # Look for H atoms bound to N
                for neighbor in atom.GetNeighbors():
                    idx = neighbor.GetIdx()
                    if neighbor.GetSymbol() == 'H' and idx in charges_dict:
                        if charges_dict[idx] > max_h_charge:
                            max_h_charge = charges_dict[idx]
                            max_h_idx = idx
    
    # Assign charges and mark the most positive H
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in charges_dict:
            atom.SetProp("charge", f"{charges_dict[idx]:.3f}")
            if is_cation and idx == max_h_idx:
                atom.SetProp("is_acidic_h", "true")

def calculate_bond_charge_differences(mol):
    """Calculate absolute charge differences for each bond"""
    bond_charge_diffs = {}
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        charge1 = float(atom1.GetProp("charge"))
        charge2 = float(atom2.GetProp("charge"))
        diff = abs(charge1 - charge2)
        bond_charge_diffs[bond.GetIdx()] = diff
    return bond_charge_diffs

def draw_molecule_with_charges(mol, bond_charge_diffs, filename="molecule_with_charges.png"):
    """Draw molecule with charge differences on bonds"""
    # Create drawing object
    drawer = Draw.MolDraw2DCairo(400, 400)
    
    # Set drawing options
    opts = drawer.drawOptions()
    opts.addAtomIndices = True
    
    # Draw the molecule
    drawer.DrawMolecule(mol)
    
    # Get the drawing
    drawer.FinishDrawing()
    
    # Save the image
    drawer.WriteDrawingText(filename)
    print(f"Image saved as {filename}")
    
    bond_partial_charge_diffs = {}
    # Print the charge differences
    print("\nBond charge differences:")
    print("(Acidic H+ bonds marked with *)")
    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        if idx in bond_charge_diffs:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # Skip regular H bonds but keep acidic H bonds
            if (atom1.GetSymbol() == 'H' and not atom1.HasProp("is_acidic_h")) or \
               (atom2.GetSymbol() == 'H' and not atom2.HasProp("is_acidic_h")):
                continue
            
            bond_charge_string = f"{atom1.GetSymbol()}{atom1.GetIdx()}-{atom2.GetSymbol()}{atom2.GetIdx()}"
            # Add asterisk if either atom is the acidic H
            if (atom1.HasProp("is_acidic_h") or atom2.HasProp("is_acidic_h")):
                bond_charge_string += "*"
            
            bond_charge_diff = bond_charge_diffs[idx]
            bond_partial_charge_diffs[bond_charge_string] = bond_charge_diff

    # Sort and print bond charge differences
    sorted_bond_partial_charge_diffs = dict(sorted(bond_partial_charge_diffs.items(), key=lambda item: item[1], reverse=True))
    for key, value in sorted_bond_partial_charge_diffs.items():
        print(f"{key}: {value:.3f}")

if __name__ == "__main__":
    # Define directory and filenames
    dirname = "data/SP_PBE0_def2TZVPP/id_108_sp_only/id108_protonated"
    charges_file = os.path.join(dirname, "molecule.out")
    xyz_file = os.path.join(dirname, "molecule.xyz")
    
    # SMILES for the molecule
    smiles = "[NH+](C(C)C)(C(C)C)CCO"
    
    # Check if molecule is cationic
    is_cation = is_cationic(dirname) or '+' in smiles
    if is_cation:
        print("\nDetected cationic/protonated species")
    
    # Read charges and XYZ data
    charges = read_charges(charges_file)
    xyz_content = read_xyz(xyz_file)
    
    if charges and xyz_content:
        # Process the molecule
        mol = xyz_to_mol(xyz_content, smiles, is_cation)
        if mol:
            assign_charges(mol, charges, is_cation)
            bond_charge_diffs = calculate_bond_charge_differences(mol)

            # Generate SMILES with explicit hydrogens
            result_smiles = Chem.MolToSmiles(mol, allHsExplicit=True)
            print(f"SMILES representation: {result_smiles}")

            # Create visualization and print bond differences
            draw_molecule_with_charges(mol, bond_charge_diffs)