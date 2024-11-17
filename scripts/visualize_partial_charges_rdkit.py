from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import rdkit

def smiles_to_mol(smiles):
    """Convert SMILES to RDKit mol object with 3D coordinates"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Could not parse SMILES")
        return None

    
    # Generate 3D coordinates
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        print("Warning: Could not generate 3D coordinates, proceeding with 2D")
        AllChem.Compute2DCoords(mol)
    
    return mol

def calculate_gasteiger_charges(mol):
    """Calculate Gasteiger partial charges for the molecule"""
    try:
        # Calculate Gasteiger charges
        rdkit.Chem.AllChem.ComputeGasteigerCharges(mol)
        
        # Store charges in atom properties
        charges = {}
        for atom in mol.GetAtoms():
            charge = float(atom.GetProp('_GasteigerCharge'))
            atom.SetProp("charge", f"{charge:.3f}")
            charges[atom.GetIdx()] = charge
            
        return charges
    except:
        print("Error: Could not calculate Gasteiger charges")
        return None

def is_cationic(smiles):
    """Check if molecule is cationic based on SMILES"""
    return '+' in smiles

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
    
    # Print the charge differences
    print("\nBond charge differences:")
    bond_partial_charge_diffs = {}
    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        if idx in bond_charge_diffs:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            
            # Skip bonds involving non-acidic hydrogens
            if (atom1.GetSymbol() == 'H' or atom2.GetSymbol() == 'H'):
                continue
            
            bond_charge_string = f"{atom1.GetSymbol()}{atom1.GetIdx()}-{atom2.GetSymbol()}{atom2.GetIdx()}"
            bond_charge_diff = bond_charge_diffs[idx]
            bond_partial_charge_diffs[bond_charge_string] = bond_charge_diff

    # Sort and print bond charge differences
    sorted_bond_partial_charge_diffs = dict(sorted(bond_partial_charge_diffs.items(), 
                                                 key=lambda item: item[1], 
                                                 reverse=True))
    for key, value in sorted_bond_partial_charge_diffs.items():
        print(f"{key}: {value:.3f}")

def analyze_molecule_charges(smiles):
    """Main function to analyze charges from SMILES"""
    print(f"\nAnalyzing molecule: {smiles}")
    
    # Check if molecule is cationic
    if is_cationic(smiles):
        print("Detected cationic species")
    
    # Convert SMILES to molecule
    mol = smiles_to_mol(smiles)
    if mol is None:
        return
    
    # Calculate charges
    charges = calculate_gasteiger_charges(mol)
    if charges is None:
        return
    
    # Calculate and visualize bond charge differences
    bond_charge_diffs = calculate_bond_charge_differences(mol)
    draw_molecule_with_charges(mol, bond_charge_diffs)
    
    # Print atom charges
    print("\nAtom charges:")
    for atom in mol.GetAtoms():
        print(f"{atom.GetSymbol()}{atom.GetIdx()}: {float(atom.GetProp('charge')):.3f}")

if __name__ == "__main__":
    # Example usage
    test_smiles = "[NH+]=C(N)NCCCC(NC(=O)CCCN)C(=O)O"
    analyze_molecule_charges(test_smiles)