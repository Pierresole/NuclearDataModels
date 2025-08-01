# Nuclear Data Models 🔬⚛️

**High-performance C++ library with Python bindings for nuclear data evaluation and R-Matrix calculations**

![Nuclear Cross Sections](https://img.shields.io/badge/Physics-Nuclear%20Data-blue?style=for-the-badge&logo=atom&logoColor=white)
![Language](https://img.shields.io/badge/C++-23-blue?logo=cplusplus)
![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Build](https://img.shields.io/badge/CMake-Eigen3-green?style=for-the-badge)

## 🎯 Showcase: Elastic Cross Section Analysis

Visualize nuclear cross sections with publication-ready plots in just 2 lines of Python:

```python
compound_system = create_compound_from_ReichMoore('Pu-239.endf')
plot_elastic_comparison(compound_system)
```

![Elastic Cross Section Example](images/cross_section_demo.png)

*Individual spin group contributions to elastic scattering in ²³⁹Pu, demonstrating resonance structures and quantum mechanical interference effects.*

## ✨ Key Features

- 🚀 **High Performance**: C++ core with Eigen3 linear algebra for fast matrix operations
- 🐍 **Python Integration**: Seamless pybind11 bindings for interactive analysis
- 📊 **Rich Visualization**: Publication-ready plots with matplotlib integration
- 🔬 **Multi-format Support**: ENDF-6 (LRF=3, LRF=7), MLBW, and R-Matrix Limited formats
- 🎛️ **Sensitivity Analysis**: Complete parameter derivative calculations for uncertainty quantification  
- ⚡ **Vectorized Operations**: Efficient energy-dependent cross section calculations
- 🧮 **Multigroup Methods**: Energy averaging for reactor physics applications
- 📈 **Resonance Analysis**: Automated peak identification and characterization
- 🔄 **Spin Group Decomposition**: Individual J-π contributions and interference effects
- 📚 **Comprehensive Documentation**: Jupyter notebooks with physics explanations

## 🏗️ Architecture Overview

Our modular design enables physicists to build complex nuclear calculations from fundamental components:

```mermaid
graph TB
    %% Main Components
    CS[CompoundSystem]:::main
    SG[SpinGroup]:::main
    CH[Channel]:::main
    PP[ParticlePair]:::main
    
    %% Supporting Classes
    RWF[RadialWaveFunctions]:::support
    RES[Resonance]:::support
    CALC[CrossSectionCalculator]:::support
    SENS[ParameterSensitivity]:::support
    
    %% Data Structures
    FP[FormalismParameters]:::data
    MLBW[MLBWParameters]:::data
    RM[RMatrixParameters]:::data
    
    %% Relationships
    CS --> SG
    SG --> CH
    SG --> RES
    CH --> PP
    CH --> RWF
    CS --> CALC
    CS --> SENS
    
    FP --> MLBW
    FP --> RM
    MLBW --> SG
    
    %% Styling
    classDef main fill:#e1f5fe,stroke:#01579b,stroke-width:3px,color:#000
    classDef support fill:#f3e5f5,stroke:#4a148c,stroke-width:2px,color:#000
    classDef data fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px,color:#000
```

### 🧮 Core Physics Classes

### 🧮 Core Physics Classes

```mermaid
classDiagram
    class CompoundSystem {
        -vector~SpinGroup~ spinGroups_
        -ParticlePair entranceParticlePair_
        +elasticCrossSection(energy) double
        +captureCrossSection(energy) double
        +fissionCrossSection(energy) double
        +totalCrossSection(energy) double
        +spinGroupElasticCrossSection(index, energy) double
        +printSpinGroupInfo() void
    }
    
    class SpinGroup {
        -double J_ "Total Angular Momentum"
        -int PJ_ "Parity"
        -vector~Channel~ channels_
        -vector~Resonance~ resonances_
        -vector~bool~ eliminatedChannels_
        +addChannel(channel) void
        +addResonance(resonance) void
        +crossSection(energy, entrancePP) double
        +elasticCrossSection(energy, entrancePP) double
        +getRetainedChannels() vector~Channel~
        +getEliminatedChannels() vector~Channel~
    }
    
    class Channel {
        -ParticlePair particlePair_
        -unsigned int l_ "Orbital Angular Momentum"
        -double effectiveChannelRadius_
        -double trueChannelRadius_
        -double spinChannel_
        +computeChannelQuantities(E, entrancePP) ChannelQuantities
        +PenetrationFactor(energy, l, entrancePP) double
        +ShiftFactor(energy, l, entrancePP) double
        +PhaseShift(energy, l, entrancePP) double
    }
    
    class ParticlePair {
        -double massParticleA_
        -double massParticleB_
        -double spinParticleA_
        -double spinParticleB_
        -double QI_ "Q-value"
        -int parityParticleA_
        -int parityParticleB_
        -int MT_ "Reaction Type"
        +k2(energy, entrancePP) double
        +reducedMass() double
        +neutronIncident(massB, spinB, QI, parityB, MT)$ ParticlePair
    }
    
    class RadialWaveFunctions {
        -double AP "Channel Radius"
        -vector~double~ P_ "Penetration Factors"
        -vector~double~ S_ "Shift Factors"
        -vector~double~ cos_phi_ "Phase Cosines"
        -vector~double~ sin_phi_ "Phase Sines"
        +P_S_Phi(l, rho2, eta) void
        +getP(l) double
        +getS(l) double
        +channelPenetrationAndShift(l, rho2) double
    }
    
    class Resonance {
        -double energy_
        -vector~double~ gamma_ "Reduced Width Amplitudes"
        +getEnergy() double
        +getGamma() vector~double~
    }
    
    class ParameterSensitivity {
        -CompoundSystem baseSystem_
        -vector~ParameterInfo~ parameters_
        -double defaultPerturbation_
        +extractParameters() void
        +computeDerivative(paramIndex, energy, crossSectionFunc) double
        +totalCrossSectionDerivatives(energy) vector~double~
        +multigroupSensitivityMatrix(energyGroups) map
        +createUniformEnergyGroups(min, max, nGroups)$ vector~EnergyGroup~
    }

    %% Relationships
    CompoundSystem *-- SpinGroup : contains
    SpinGroup *-- Channel : contains
    SpinGroup *-- Resonance : contains
    Channel *-- ParticlePair : uses
    Channel ..> RadialWaveFunctions : computes with
    CompoundSystem ..> ParameterSensitivity : analyzed by
    
    %% Nested structs
    class ChannelQuantities {
        +double P "Penetration"
        +double S "Shift"
        +double phi "Phase Shift"
        +double k2 "Wave Number²"
    }
    
    class ParameterInfo {
        +Type type
        +size_t spinGroupIndex
        +size_t resonanceIndex
        +string description
        +double nominalValue
    }
    
    Channel --> ChannelQuantities : returns
    ParameterSensitivity --> ParameterInfo : uses
```

## 🚀 Quick Start

### Prerequisites
- **C++23** compatible compiler
- **CMake** 3.15+
- **Python** 3.8+ with development headers
- **Eigen3** (automatically fetched)
- **pybind11** (automatically fetched)

### Build & Install

```bash
# Clone and setup
git clone https://github.com/YourUsername/NuclearDataModels.git
cd NuclearDataModels

# Build the library
mkdir build && cd build
cmake -D Python3_EXECUTABLE=$(which python3) ..
make -j$(nproc)

# Test the installation
cd ../notebooks/RMatrix
python -c "import sys; sys.path.append('../../build/python'); import pyRMatrix; print('✅ Import successful!')"
```

### First Calculation

```python
import sys
sys.path.append('build/python')
import pyRMatrix
from scripts.compoundFromENDFtk import create_compound_from_ReichMoore

# Load nuclear data (Pu-239 example)
compound = create_compound_from_ReichMoore('n-094_Pu_239.endf')

# Calculate cross sections
energy = 1.0  # eV
print(f"Elastic:  {compound.elasticCrossSection(energy):.2f} barns")
print(f"Capture:  {compound.captureCrossSection(energy):.2f} barns") 
print(f"Fission:  {compound.fissionCrossSection(energy):.2f} barns")
print(f"Total:    {compound.totalCrossSection(energy):.2f} barns")

# Quick visualization
from scripts.plot_cross_sections import plot_all_reactions
plot_all_reactions(compound)
```

## 📊 Advanced Features

### 🎛️ Parameter Sensitivity Analysis

Compute numerical derivatives of cross sections with respect to all resonance parameters:

```python
# Initialize sensitivity analyzer
sensitivity = pyRMatrix.ParameterSensitivity(compound, perturbation=1e-6)

# Analyze sensitivity at thermal energy
energy = 0.0253  # eV (thermal)
derivatives = sensitivity.totalCrossSectionDerivatives(energy)

# Multigroup sensitivity matrix
energy_groups = sensitivity.createLogarithmicEnergyGroups(0.01, 100, 10)
sensitivity_matrix = sensitivity.computeMultigroupSensitivityMatrix(energy_groups)

print(f"Computed {len(derivatives)} parameter sensitivities")
print(f"Most sensitive parameter: {max(derivatives):.2e}")
```

### 🔬 Interactive Jupyter Notebooks

Explore the full capabilities with our comprehensive notebook examples:

#### [`resonan_from_endfTK_LRF7.ipynb`](notebooks/RMatrix/resonan_from_endfTK_LRF7.ipynb) - Main Demonstration
- 📈 **Cross Section Visualization**: Publication-ready plots of elastic, capture, and fission cross sections
- 🔍 **Spin Group Analysis**: Individual J-π group contributions and interference effects  
- 📊 **Parameter Sensitivity**: Numerical derivatives and uncertainty quantification
- ⚡ **Performance Comparison**: Benchmark different computational approaches
- 🎯 **Resonance Identification**: Automated peak finding and resonance characterization

#### [`sensitivity.ipynb`](notebooks/RMatrix/sensitivity.ipynb) - Advanced Analysis
- 🧮 **Multigroup Calculations**: Energy-averaged cross sections for reactor physics
- 📉 **Sensitivity Matrices**: Complete parameter correlation analysis
- 🎚️ **Perturbation Studies**: Understanding parameter impact on observables

#### Key Notebook Features:
```python
# One-liner nuclear data analysis
compound_system = create_compound_from_ReichMoore('Pu-239.endf')
plot_elastic_comparison(compound_system)  # → Beautiful publication plot!

# Advanced spin group breakdown  
plot_spin_group_breakdown(compound_system, 'fission')

# Parameter sensitivity heatmaps
sensitivity.computeMultigroupSensitivityMatrix(energy_groups)
```

### 🔄 Multi-format Input Support

```python
# Reich-Moore format (LRF=3)
compound_rm = create_compound_from_ReichMoore('plutonium.endf')

# R-Matrix Limited format (LRF=7) 
compound_rml = create_compound_from_RMatrix('copper.endf')

# MLBW format support
wrapper = pyRMatrix.FormalismParametersWrapper()
wrapper.fill(endf_data_dict)
mlbw = pyRMatrix.MLBW(wrapper.data.parameters.mlbw)
```

### 📈 Publication-Ready Plotting

```python
from scripts.plot_cross_sections import *

# Individual spin group contributions
plot_spin_group_breakdown(compound, 'elastic')

# Compare all reaction types
plot_all_reactions(compound)

# Customizable energy ranges and styling
plot_elastic_comparison(compound, energy_range=(-3, 3), save_path='figure.png')
```

## 🧪 Example Calculations

### Resonance Analysis
```python
# Examine individual resonances
compound.printSpinGroupInfo()

# Calculate at resonance energy
resonance_energy = 0.296  # eV (Pu-239 resonance)
xs_at_resonance = compound.totalCrossSection(resonance_energy)
print(f"Cross section at resonance: {xs_at_resonance:.0f} barns")
```

### Energy-dependent Cross Sections
```python
import numpy as np

# Thermal to fast neutron range
energies = np.logspace(-3, 6, 1000)  # 1 meV to 1 MeV
cross_sections = [compound.totalCrossSection(E) for E in energies]

# Find resonance peaks
peaks = np.where(np.array(cross_sections) > 1000)[0]
print(f"Found {len(peaks)} resonance peaks above 1000 barns")
```

## 🎓 Physical Background

This library implements modern nuclear data evaluation methods based on:

- **R-Matrix Theory**: Describes nuclear reactions through scattering matrix formalism
- **Multi-Level Breit-Wigner (MLBW)**: Resonance parameterization for isolated resonances  
- **Reich-Moore**: Advanced formalism handling overlapping resonances
- **Radial Wave Functions**: Quantum mechanical boundary conditions at nuclear surface

### Key References
- 📖 **SAMMY Manual** - Oak Ridge National Laboratory
- 📖 **G. Ferran's Thesis** - Advanced R-Matrix implementations
- 📖 **JEFF Report 18** - European nuclear data evaluation guidelines
- 📖 **ENDF Manual** - International nuclear data format standards
- 📖 **Neutron Interaction Theory** - A. Foderaro

## 🏗️ Project Structure

```
NuclearDataModels/
├── 📁 include/RMatrix/          # C++ header files
│   ├── CompoundSystem.h         # Main calculation engine
│   ├── SpinGroup.h             # Quantum number groupings
│   ├── Channel.h               # Reaction channels
│   ├── ParticlePair.h          # Particle definitions
│   ├── RadialWaveFunctions.h   # Quantum wave functions
│   └── ParameterSensitivity.h  # Uncertainty quantification
├── 📁 src/RMatrix/             # C++ implementation
├── 📁 bindings/python/         # Pybind11 interfaces  
├── 📁 notebooks/RMatrix/       # Jupyter examples
│   ├── resonan_from_endfTK_LRF7.ipynb  # Main demo
│   ├── resonance_from_tk.ipynb         # Import examples
│   └── sensitivity.ipynb               # Sensitivity analysis
├── 📁 scripts/                 # Python utilities
│   ├── compoundFromENDFtk.py   # ENDF file parsers
│   └── plot_cross_sections.py # Visualization tools
└── 📁 build/                   # CMake build output
    └── python/pyRMatrix.so     # Compiled Python module
```

## 🤝 Contributing

We welcome contributions from the nuclear physics community! Please see our [Contributing Guidelines](CONTRIBUTING.md) for:

- 🐛 Bug reports and feature requests
- 🔬 Physics model improvements  
- 📊 Visualization enhancements
- 📖 Documentation and examples
- 🧪 Test cases and benchmarks

## 📄 License

This project is licensed under the [LICENSE](LICENSE) - see the file for details.

## 🙏 Acknowledgments

Built upon decades of nuclear physics research and computational methods development. Special thanks to the nuclear data evaluation communities at ORNL, CEA, and IAEA for their foundational work in R-Matrix theory and implementation.

---

<div align="center">

**⚛️ Advancing Nuclear Science Through Computational Excellence ⚛️**

[![GitHub](https://img.shields.io/badge/GitHub-View%20Source-black?style=for-the-badge&logo=github)](https://github.com/YourUsername/NuclearDataModels)
[![Documentation](https://img.shields.io/badge/Docs-Read%20More-blue?style=for-the-badge&logo=readthedocs)](https://your-docs-url.com)

</div>
