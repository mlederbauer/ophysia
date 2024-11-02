from enum import Enum
from typing import Dict, List

class Functional(Enum):
    B3LYP = "B3LYP"
    PBE0 = "PBE0"
    wB97X = "wB97X"

class FunctionalManager:
    @staticmethod
    def get_functional_keywords(functional: Functional) -> List[str]:
        """Get ORCA keywords for functional.

        Args:
            functional (Functional): Functional type
        
        Returns:
            List[str]: List of ORCA keywords
        
        Raises:
            KeyError: If functional type not found
        """

        if not isinstance(functional, Functional):
            raise ValueError(f"Unknown functional: {functional}")

        keywords = {
            Functional.B3LYP: ["B3LYP"],
            Functional.PBE0: ["PBE0"],
            Functional.wB97X: ["wB97X"],
        }
        return keywords.get(functional, [])