from enum import Enum
from typing import Dict, List, Optional
from dataclasses import dataclass

@dataclass
class BasisSetInfo:
    """Information about a basis set."""
    keyword: str
    description: str
    elements: Optional[List[str]] = None  # List of supported elements, None means all
    has_ecp: bool = False  # Whether basis set includes effective core potential
    is_auxiliary: bool = False  # Whether it's an auxiliary basis for RI/RIJCOSX

class BasisSet(Enum):
    """Common basis sets available in ORCA."""
    # Minimal basis sets
    STO_3G = "STO-3G"
    MINI = "MINI"
    
    # Split valence basis sets
    def2_SV = "def2-SV"
    def2_SVP = "def2-SVP"
    def2_SVPD = "def2-SVPD"
    
    # Triple zeta basis sets
    def2_TZVP = "def2-TZVP"
    def2_TZVPP = "def2-TZVPP"
    def2_TZVPD = "def2-TZVPD"
    
    # Quadruple zeta basis sets
    def2_QZVP = "def2-QZVP"
    def2_QZVPP = "def2-QZVPP"
    
    # Dunning correlation consistent basis sets
    cc_pVDZ = "cc-pVDZ"
    cc_pVTZ = "cc-pVTZ"
    cc_pVQZ = "cc-pVQZ"
    aug_cc_pVDZ = "aug-cc-pVDZ"
    aug_cc_pVTZ = "aug-cc-pVTZ"
    
    # Special purpose basis sets
    def2_UNIVERSAL = "def2-UNIVERSAL"
    SARC_DKH = "SARC-DKH"

class BasisSetManager:
    """Manages basis set information and keywords."""
    
    _BASIS_INFO: Dict[BasisSet, BasisSetInfo] = {
        BasisSet.STO_3G: BasisSetInfo(
            keyword="STO-3G",
            description="Minimal basis set, suitable for quick preliminary calculations"
        ),
        BasisSet.def2_SVP: BasisSetInfo(
            keyword="def2-SVP",
            description="Split valence polarized basis set, good for initial geometry optimizations"
        ),
        BasisSet.def2_TZVP: BasisSetInfo(
            keyword="def2-TZVP",
            description="Triple zeta valence polarized, good quality for production calculations"
        ),
        BasisSet.def2_QZVP: BasisSetInfo(
            keyword="def2-QZVP",
            description="Quadruple zeta valence polarized, high accuracy calculations"
        ),
        BasisSet.cc_pVDZ: BasisSetInfo(
            keyword="cc-pVDZ",
            description="Correlation consistent double zeta basis set"
        ),
        BasisSet.aug_cc_pVTZ: BasisSetInfo(
            keyword="aug-cc-pVTZ",
            description="Augmented triple zeta basis set, good for anions and excited states"
        ),
    }

    @classmethod
    def get_basis_keywords(cls, basis: BasisSet) -> List[str]:
        """Get ORCA keywords for a basis set.

        Args:
            basis (BasisSet): Basis set
        
        Returns:
            List[str]: List of ORCA keywords
        
        Raises:
            ValueError: If basis set not found
        """

        try:
            info = cls._BASIS_INFO[basis]
            return [f"! {info.keyword}"]
        except KeyError:
            raise ValueError(f"Unknown basis set: {basis}")

    @classmethod
    def get_basis_info(cls, basis: BasisSet) -> BasisSetInfo:
        """Get detailed information about a basis set.

        Args:
            basis (BasisSet): Basis set
        
        Returns:
            BasisSetInfo: Information about the basis set
        
        Raises:
            ValueError: If basis set not found
        """

        try:
            return cls._BASIS_INFO[basis]
        except KeyError:
            raise ValueError(f"Unknown basis set: {basis}")

    @classmethod
    def get_auxiliary_basis(cls, basis: BasisSet) -> Optional[str]:
        """Get appropriate auxiliary basis set for RI/RIJCOSX calculations.

        Args:
            basis (BasisSet): Main basis set

        Returns:
            Optional[str]: Auxiliary basis set keyword, or None if not needed
        """

        auxiliary_map = {
            BasisSet.def2_SVP: "/J def2/J",
            BasisSet.def2_TZVP: "/J def2/J",
            BasisSet.def2_QZVP: "/J def2/J",
        }
        return auxiliary_map.get(basis)

    @classmethod
    def suggest_basis(cls, 
                     elements: List[str], 
                     accuracy: str = "medium",
                     diffuse: bool = False) -> BasisSet:
        """Suggest appropriate basis set based on elements and desired accuracy.

        Args:
            elements (List[str]): List of atomic symbols
            accuracy (str): Desired accuracy level (low, medium, high)
            diffuse (bool): Whether diffuse functions are needed
        
        Returns:
            BasisSet: Recommended basis set
        """

        accuracy_map = {
            "low": BasisSet.def2_SVP,
            "medium": BasisSet.def2_TZVP,
            "high": BasisSet.def2_QZVP
        }
        
        # Check if any heavy elements present (Z > 36)
        heavy_elements = {"Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", 
                         "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"}
        
        if any(elem in heavy_elements for elem in elements):
            return BasisSet.SARC_DKH
        
        suggested = accuracy_map.get(accuracy.lower(), BasisSet.def2_TZVP)
        
        # If diffuse functions needed, switch to appropriate basis
        if diffuse:
            if suggested == BasisSet.def2_SVP:
                return BasisSet.def2_SVPD
            elif suggested == BasisSet.def2_TZVP:
                return BasisSet.def2_TZVPD
            
        return suggested
