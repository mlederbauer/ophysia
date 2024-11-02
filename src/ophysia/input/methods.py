from enum import Enum
from typing import List

class CalculationType(Enum):
    SINGLE_POINT = "SP"
    SINGLE_POINT_CUSTOM = "SP-CUSTOM"
    OPT = "OPT"
    FREQ = "FREQ"
    OPT_FREQ_GOAT = "OPT-FREQ-GOAT"
    OPT_CUSTOM = "OPT-CUSTOM"
    D3 = "D3"
    D4 = "D4"

class MethodManager:
    @staticmethod
    def get_calculation_keywords(calc_type: CalculationType) -> List[str]:
        """Get ORCA keywords for calculation type.

        Args:
            calc_type (CalculationType): Calculation type
        
        Returns:
            List[str]: List of ORCA keywords
    
        Raises:
            KeyError: If calculation type not found
        """

        if not isinstance(calc_type, CalculationType):
            raise ValueError(f"Unknown calculation type: {calc_type}")

        keywords = {
            CalculationType.SINGLE_POINT: ["SP"],
            CalculationType.SINGLE_POINT_CUSTOM: ["SP RIJCOSX def2/J TIGHTSCF"],
            CalculationType.OPT: ["OPT"],
            CalculationType.FREQ: ["FREQ"],
            CalculationType.OPT_FREQ_GOAT: ["OPT FREQ GOAT"],
            CalculationType.OPT_CUSTOM: ["OPT RIJCOSX FREQ def2/J TIGHTOPT"],
        }
        return keywords.get(calc_type, [])
    
    @staticmethod
    def get_dispersion_keywords(calc_type: CalculationType) -> List[str]:
        """Get ORCA keywords for dispersion correction.

        Args:
            calc_type (CalculationType): Calculation type
        
        Returns:
            List[str]: List of ORCA keywords
    
        Raises:
            KeyError: If calculation type not found
        """

        if not isinstance(calc_type, CalculationType):
            raise ValueError(f"Unknown calculation type: {calc_type}")

        keywords = {
            CalculationType.D3: ["D3"],
            CalculationType.D4: ["D4"],
        }
        return keywords.get(calc_type, [])