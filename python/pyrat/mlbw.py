# mlbw.py
from .resonance_data import ResonanceDataMLBW
import pyrat_bindings  # Import the bindings module for C++ interaction

class MLBW:
    def __init__(self, data: ResonanceDataMLBW):
        self.data = data
        # Initialize the C++ MLBW object with data if needed
        self.cpp_obj = pyrat_bindings.MLBW(self.data)

    def recon(self):
        # Call the C++ recon function
        return self.cpp_obj.recon()
