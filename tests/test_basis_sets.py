import pytest
from ophysia.input.basis_sets import BasisSet, BasisSetManager, BasisSetInfo

def test_basis_set_enum_values():
    """Test that basis set enum contains expected values."""
    assert BasisSet.def2_SVP.value == "def2-SVP"
    assert BasisSet.def2_TZVP.value == "def2-TZVP"
    assert BasisSet.cc_pVDZ.value == "cc-pVDZ"

def test_get_basis_keywords():
    """Test basis set keyword generation."""
    manager = BasisSetManager()
    
    # Test def2-SVP keywords
    svp_keywords = manager.get_basis_keywords(BasisSet.def2_SVP)
    assert isinstance(svp_keywords, list)
    assert "! def2-SVP" in svp_keywords
    
    # Test def2-TZVP keywords
    tzvp_keywords = manager.get_basis_keywords(BasisSet.def2_TZVP)
    assert isinstance(tzvp_keywords, list)
    assert "! def2-TZVP" in tzvp_keywords

def test_get_basis_info():
    """Test retrieval of basis set information."""
    manager = BasisSetManager()
    info = manager.get_basis_info(BasisSet.def2_SVP)
    
    assert isinstance(info, BasisSetInfo)
    assert "split valence" in info.description.lower()

def test_invalid_basis():
    """Test handling of invalid basis set."""
    manager = BasisSetManager()
    with pytest.raises(ValueError):
        manager.get_basis_keywords(BasisSet("INVALID"))

def test_get_auxiliary_basis():
    """Test auxiliary basis set retrieval."""
    manager = BasisSetManager()
    aux_basis = manager.get_auxiliary_basis(BasisSet.def2_SVP)
    assert aux_basis == "/J def2/J"

@pytest.mark.parametrize("elements,accuracy,diffuse,expected", [
    (["H", "C", "N", "O"], "low", False, BasisSet.def2_SVP),
    (["H", "C", "N", "O"], "medium", False, BasisSet.def2_TZVP),
    (["H", "C", "N", "O"], "high", False, BasisSet.def2_QZVP),
    (["H", "C", "N", "O"], "low", True, BasisSet.def2_SVPD),
    (["Ag", "H"], "medium", False, BasisSet.SARC_DKH),
])
def test_suggest_basis(elements, accuracy, diffuse, expected):
    """Test basis set suggestion system."""
    manager = BasisSetManager()
    suggested = manager.suggest_basis(elements, accuracy, diffuse)
    assert suggested == expected