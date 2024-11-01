import pytest
from rdkit import Chem
from pathlib import Path
from ophysia.input import functionals, methods, geometry, utils, basis_sets

# Test data
EXAMPLE_XYZ = """3
Water molecule
O  0.0000  0.0000  0.0000
H  0.0000  0.7580  0.5040
H  0.0000  -0.7580  0.5040
"""

@pytest.fixture
def xyz_file(tmp_path):
    """Create a temporary XYZ file for testing."""
    xyz_path = tmp_path / "water.xyz"
    xyz_path.write_text(EXAMPLE_XYZ)
    return xyz_path

# Test Functionals
def test_functional_enum_values():
    """Test that functional enum contains expected values."""
    assert functionals.Functional.B3LYP.value == "B3LYP"
    assert functionals.Functional.PBE0.value == "PBE0"
    assert functionals.Functional.wB97X.value == "wB97X"

def test_get_functional_keywords():
    """Test functional keyword generation."""
    manager = functionals.FunctionalManager()
    
    # Test B3LYP keywords
    b3lyp_keywords = manager.get_functional_keywords(functionals.Functional.B3LYP)
    assert isinstance(b3lyp_keywords, list)
    assert "! B3LYP" in b3lyp_keywords
    
    # Test PBE0 keywords
    pbe0_keywords = manager.get_functional_keywords(functionals.Functional.PBE0)
    assert isinstance(pbe0_keywords, list)
    assert "! PBE0" in pbe0_keywords

def test_invalid_functional():
    """Test handling of invalid functional."""
    manager = functionals.FunctionalManager()
    with pytest.raises(ValueError):
        manager.get_functional_keywords(functionals.Functional("INVALID"))

# Test Methods
def test_calculation_type_enum_values():
    """Test that calculation type enum contains expected values."""
    assert methods.CalculationType.SINGLE_POINT.value == "SP"
    assert methods.CalculationType.OPT.value == "OPT"
    assert methods.CalculationType.FREQ.value == "FREQ"

def test_get_calculation_keywords():
    """Test calculation keyword generation."""
    manager = methods.MethodManager()
    
    # Test single point keywords
    sp_keywords = manager.get_calculation_keywords(methods.CalculationType.SINGLE_POINT)
    assert isinstance(sp_keywords, list)
    assert "! SP" in sp_keywords
    
    # Test optimization keywords
    opt_keywords = manager.get_calculation_keywords(methods.CalculationType.OPT)
    assert isinstance(opt_keywords, list)
    assert "! OPT" in opt_keywords

def test_invalid_calculation_type():
    """Test handling of invalid calculation type."""
    manager = methods.MethodManager()
    with pytest.raises(ValueError):
        manager.get_calculation_keywords(methods.CalculationType("INVALID"))

# Test Geometry
def test_geometry_builder_initialization(xyz_file):
    """Test GeometryBuilder initialization."""
    builder = geometry.GeometryBuilder(xyz_file)
    assert isinstance(builder.xyz_file, Path)
    assert builder.xyz_file.exists()

def test_create_input_geometry(xyz_file):
    """Test creation of ORCA geometry input."""
    builder = geometry.GeometryBuilder(xyz_file)
    geom_lines = builder.create_input_geometry(charge=0, multiplicity=1)
    
    assert isinstance(geom_lines, list)
    assert len(geom_lines) > 0
    assert geom_lines[0] == "* xyz 0 1"
    assert "O" in geom_lines[1]  # Check atom symbol
    assert "*" in geom_lines[-1]  # Check closing marker

def test_geometry_with_charge_multiplicity(xyz_file):
    """Test geometry creation with different charge and multiplicity."""
    builder = geometry.GeometryBuilder(xyz_file)
    geom_lines = builder.create_input_geometry(charge=1, multiplicity=2)
    
    assert geom_lines[0] == "* xyz 1 2"

def test_invalid_xyz_file():
    """Test handling of invalid XYZ file."""
    with pytest.raises(FileNotFoundError):
        builder = geometry.GeometryBuilder(Path("nonexistent.xyz"))
        builder.create_input_geometry()

# Test Utils
@pytest.mark.skipif(not utils.is_obabel_available(), 
                   reason="OpenBabel not available")
def test_smiles_to_xyz(tmp_path):
    """Test SMILES to XYZ conversion."""
    output_file = tmp_path / "molecule.xyz"
    smiles = "O"  # Water
    
    result = utils.smiles_to_xyz(smiles, output_file)
    
    assert result.exists()
    assert result.suffix == ".xyz"
    
    # Check content
    content = result.read_text()
    assert "O" in content  # Should contain oxygen atom

@pytest.mark.skipif(not utils.is_obabel_available(), 
                   reason="OpenBabel not available")
def test_invalid_smiles(tmp_path):
    """Test handling of invalid SMILES."""
    output_file = tmp_path / "invalid.xyz"
    smiles = "C1CCO"
    
    with pytest.raises(ValueError):
        utils.smiles_to_xyz(smiles, output_file)


def test_stereochemistry_handling(tmp_path):
    """Test handling of stereochemistry in SMILES."""
    output_file = tmp_path / "chiral.xyz"
    chiral_smiles = "C[C@H](O)C(=O)O"  # Lactic acid
    
    result = utils.smiles_to_xyz(chiral_smiles, output_file)
    assert result.exists()
    
    # Read the content
    content = result.read_text().split('\n')
    natoms = int(content[0])
    
    # Count actual coordinate lines (excluding empty lines and comments)
    coord_lines = [line.strip() for line in content[2:] if line.strip()]
    assert len(coord_lines) == natoms
    
    # Verify that the number of coordinate lines matches the expected number
    mol = Chem.MolFromSmiles(chiral_smiles)
    # Note: OpenBabel adds explicit hydrogens while RDKit's GetAtoms() only counts heavy atoms
    # We should compare the coordinate lines count with the total atom count including hydrogens
    total_atoms = len(list(mol.GetAtoms())) + sum(atom.GetTotalNumHs() for atom in mol.GetAtoms()) # type: ignore[no-untyped-call, call-arg]
    assert natoms == total_atoms

# Integration Tests
def test_full_input_generation(xyz_file, tmp_path):
    """Test full input file generation process."""
    # Setup
    builder = geometry.GeometryBuilder(xyz_file)
    func_manager = functionals.FunctionalManager()
    method_manager = methods.MethodManager()
    
    # Get components
    geom_lines = builder.create_input_geometry()
    func_keywords = func_manager.get_functional_keywords(functionals.Functional.B3LYP)
    method_keywords = method_manager.get_calculation_keywords(methods.CalculationType.OPT)
    
    # Combine components
    input_lines = [
        "# ORCA Input File",
        *method_keywords,
        *func_keywords,
        "! def2-SVP",
        "",
        *geom_lines
    ]
    
    # Write input file
    input_file = tmp_path / "test.inp"
    with open(input_file, "w") as f:
        f.write("\n".join(input_lines))
    
    # Verify
    assert input_file.exists()
    content = input_file.read_text()
    assert "! OPT" in content
    assert "! B3LYP" in content
    assert "* xyz" in content

# Parameterized Tests
@pytest.mark.parametrize("functional,expected_keyword", [
    (functionals.Functional.B3LYP, "! B3LYP"),
    (functionals.Functional.PBE0, "! PBE0"),
    (functionals.Functional.wB97X, "! wB97X"),
])
def test_functional_keywords_parametrized(functional, expected_keyword):
    """Test various functional keywords."""
    manager = functionals.FunctionalManager()
    keywords = manager.get_functional_keywords(functional)
    assert expected_keyword in keywords

@pytest.mark.parametrize("calc_type,expected_keyword", [
    (methods.CalculationType.SINGLE_POINT, "! SP"),
    (methods.CalculationType.OPT, "! OPT"),
    (methods.CalculationType.FREQ, "! FREQ"),
])
def test_calculation_keywords_parametrized(calc_type, expected_keyword):
    """Test various calculation keywords."""
    manager = methods.MethodManager()
    keywords = manager.get_calculation_keywords(calc_type)
    assert expected_keyword in keywords

# Configuration Tests
def test_basis_sets():
    """Test basis set handling."""
    # Implement once basis_sets.py is created
    pass

def test_method_combinations():
    """Test combining different methods."""
    # Test combining OPT and FREQ, etc.
    pass

# Helper function for utils.py
def test_is_obabel_available():
    """Test OpenBabel availability check."""
    result = utils.is_obabel_available()
    assert isinstance(result, bool)

# Test Basis Sets
def test_basis_set_enum_values():
    """Test that basis set enum contains expected values."""
    assert basis_sets.BasisSet.def2_SVP.value == "def2-SVP"
    assert basis_sets.BasisSet.def2_TZVP.value == "def2-TZVP"
    assert basis_sets.BasisSet.cc_pVDZ.value == "cc-pVDZ"

def test_get_basis_keywords():
    """Test basis set keyword generation."""
    manager = basis_sets.BasisSetManager()
    
    # Test def2-SVP keywords
    svp_keywords = manager.get_basis_keywords(basis_sets.BasisSet.def2_SVP)
    assert isinstance(svp_keywords, list)
    assert "! def2-SVP" in svp_keywords
    
    # Test def2-TZVP keywords
    tzvp_keywords = manager.get_basis_keywords(basis_sets.BasisSet.def2_TZVP)
    assert isinstance(tzvp_keywords, list)
    assert "! def2-TZVP" in tzvp_keywords

def test_get_basis_info():
    """Test retrieval of basis set information."""
    manager = basis_sets.BasisSetManager()
    info = manager.get_basis_info(basis_sets.BasisSet.def2_SVP)
    
    assert isinstance(info, basis_sets.BasisSetInfo)
    assert "split valence" in info.description.lower()

def test_invalid_basis():
    """Test handling of invalid basis set."""
    manager = basis_sets.BasisSetManager()
    with pytest.raises(ValueError):
        manager.get_basis_keywords(basis_sets.BasisSet("INVALID"))

def test_get_auxiliary_basis():
    """Test auxiliary basis set retrieval."""
    manager = basis_sets.BasisSetManager()
    aux_basis = manager.get_auxiliary_basis(basis_sets.BasisSet.def2_SVP)
    assert aux_basis == "/J def2/J"

@pytest.mark.parametrize("elements,accuracy,diffuse,expected", [
    (["H", "C", "N", "O"], "low", False, basis_sets.BasisSet.def2_SVP),
    (["H", "C", "N", "O"], "medium", False, basis_sets.BasisSet.def2_TZVP),
    (["H", "C", "N", "O"], "high", False, basis_sets.BasisSet.def2_QZVP),
    (["H", "C", "N", "O"], "low", True, basis_sets.BasisSet.def2_SVPD),
    (["Ag", "H"], "medium", False, basis_sets.BasisSet.SARC_DKH),
])
def test_suggest_basis(elements, accuracy, diffuse, expected):
    """Test basis set suggestion system."""
    manager = basis_sets.BasisSetManager()
    suggested = manager.suggest_basis(elements, accuracy, diffuse)
    assert suggested == expected