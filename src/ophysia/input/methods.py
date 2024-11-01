from enum import Enum
from typing import List

class CalculationType(Enum):
    SINGLE_POINT = "SP"
    OPT = "OPT"
    FREQ = "FREQ"

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

        keywords = {
            CalculationType.SINGLE_POINT: ["! SP"],
            CalculationType.OPT: ["! OPT"],
            CalculationType.FREQ: ["! FREQ"],
        }
        return keywords.get(calc_type, [])