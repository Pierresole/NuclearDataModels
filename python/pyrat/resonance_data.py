# resonance_data.py
from dataclasses import dataclass, field
from typing import List

@dataclass
class RPlist:
    ER: List[float] = field(default_factory=list)
    AJ: List[float] = field(default_factory=list)
    GN: List[float] = field(default_factory=list)
    GG: List[float] = field(default_factory=list)
    GFA: List[float] = field(default_factory=list)
    GFB: List[float] = field(default_factory=list)

@dataclass
class Llist:
    AWRI: float
    APL: float
    L: int
    RPlist: RPlist

@dataclass
class ResonanceDataMLBW:
    SPI: float
    AP: float
    LAD: int
    Llist: List['Llist'] = field(default_factory=list)
