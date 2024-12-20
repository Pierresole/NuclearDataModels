from typing import Any
from .mlbw import MLBW  # Import other formalism classes similarly
from .resonance_data import ResonanceDataMLBW

class Formalism:
    @staticmethod
    def extract_resonance_data(resorange: Any, formalism_type: str) -> Any:
        if formalism_type.lower() == 'mlbw':
            data = Formalism._extract_mlbw_data(resorange)
            return MLBW(data)
        # Handle other formalisms like SLBW, RM, RMatrix here
        else:
            raise ValueError(f"Unsupported formalism type: {formalism_type}")

    @staticmethod
    def _extract_mlbw_data(resorange) -> ResonanceDataMLBW:
        SPI = resorange.parameters.SPI
        AP = resorange.parameters.AP
        LAD = resorange.parameters.LAD

        Llist = []
        for i in range(resorange.parameters.NLS):
            l_value = resorange.parameters.l_values[i]

            AWRI = l_value.AWRI
            APL = l_value.APL
            L = l_value.L

            ER = [l_value.ER[j] for j in range(l_value.NRS)]
            AJ = [l_value.AJ[j] for j in range(l_value.NRS)]
            GN = [l_value.GN[j] for j in range(l_value.NRS)]
            GG = [l_value.GG[j] for j in range(l_value.NRS)]
            GFA = [l_value.GFA[j] for j in range(l_value.NRS)]
            GFB = [l_value.GFB[j] for j in range(l_value.NRS)]

            rp_list = RPlist(ER=ER, AJ=AJ, GN=GN, GG=GG, GFA=GFA, GFB=GFB)
            l_list_item = Llist(AWRI=AWRI, APL=APL, L=L, RPlist=rp_list)
            Llist.append(l_list_item)

        return ResonanceDataMLBW(SPI=SPI, AP=AP, LAD=LAD, Llist=Llist)
