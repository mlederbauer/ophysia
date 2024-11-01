from enum import Enum
from typing import List

class CalculationType(Enum):
    SINGLE_POINT = "SP"
    OPT = "OPT"
    FREQ = "FREQ"
    OPT_FREQ_GOAT = "OPT-FREQ-GOAT"
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
            CalculationType.OPT: ["OPT"],
            CalculationType.FREQ: ["FREQ"],
            CalculationType.OPT_FREQ_GOAT: ["OPT FREQ GOAT"],
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