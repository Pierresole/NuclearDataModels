# Nuclear Data Models 🔬⚛️

**High-performance C++ library with Python bindings for nuclear data evaluation and R-Matrix calculations**

![Nuclear Cross Sections](https://img.shields.io/badge/Physics-Nuclear%20Data-blue?style=for-the-badge&logo=atom&logoColor=white)
![Language](https://img.shields.io/badge/C++-23-blue?logo=cplusplus)
![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Build](https://img.shields.io/badge/CMake-Eigen3-green?style=for-the-badge)

## 🎯 Showcase: Intuitive Nuclear Physics API

Experience the elegance of our C++ bindings - nuclear physics calculations become as simple as calling methods:

```python
import numpy as np
import matplotlib.pyplot as plt
from scripts.compoundFromENDFtk import create_compound_from_ReichMoore

# Load nuclear data with one line
compound_system = create_compound_from_ReichMoore('Pu-239.endf')

# Inspect the nuclear structure
compound_system.printSpinGroupInfo()
# → Spin Group 0/2 (0.5, 1): has 3 channels.
# → Spin Group 1/2 (1.5, 1): has 3 channels.

# Calculate cross sections at any energy
thermal_energy = 0.0253  # eV
print(f"At thermal energy ({thermal_energy} eV):")
print(f"  Elastic:  {compound_system.elasticCrossSection(thermal_energy):.1f} barns")
print(f"  Capture:  {compound_system.captureCrossSection(thermal_energy):.1f} barns") 
print(f"  Fission:  {compound_system.fissionCrossSection(thermal_energy):.1f} barns")
print(f"  Total:    {compound_system.totalCrossSection(thermal_energy):.1f} barns")
```

![Nuclear Reactions](images/nuclear_reactions_demo.png)

*Multiple nuclear reaction types calculated with simple, intuitive method calls on the compound system object.*

```python
# Individual spin group contributions (J-π physics!)
j_half_elastic = compound_system.spinGroupElasticCrossSection(0, thermal_energy)
j_3half_elastic = compound_system.spinGroupElasticCrossSection(1, thermal_energy)
print(f"J=1/2 elastic: {j_half_elastic:.1f} barns")
print(f"J=3/2 elastic: {j_3half_elastic:.1f} barns")

# Energy-dependent analysis
energies = np.logspace(-2, 2, 1000)  # 0.01 to 100 eV
elastic_xs = [compound_system.elasticCrossSection(E) for E in energies]
capture_xs = [compound_system.captureCrossSection(E) for E in energies]

# Create publication-ready plot
plt.figure(figsize=(10, 6))
plt.loglog(energies, elastic_xs, label='Elastic', linewidth=2.5)
plt.loglog(energies, capture_xs, label='Capture (n,γ)', linewidth=2.5)
plt.xlabel('Energy (eV)')
plt.ylabel('Cross Section (barns)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

![Spin Groups](images/spin_groups_demo.png)

*Quantum mechanical spin group decomposition showing individual J-π contributions to elastic scattering.*

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

# Nuclear physics at your fingertips!
energy = 1.0  # eV
print(f"Elastic:  {compound.elasticCrossSection(energy):.2f} barns")
print(f"Capture:  {compound.captureCrossSection(energy):.2f} barns") 
print(f"Fission:  {compound.fissionCrossSection(energy):.2f} barns")
print(f"Total:    {compound.totalCrossSection(energy):.2f} barns")

# Spin group decomposition (quantum mechanics made easy!)
n_groups = len([sg for sg in compound.spinGroups()])
for i in range(n_groups):
    J = compound.getSpinGroup(i).getJ()
    elastic_j = compound.spinGroupElasticCrossSection(i, energy)
    print(f"J={J} elastic: {elastic_j:.2f} barns")

# Quick visualization with our plotting utilities
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

While the compound system API provides direct access to all physics calculations, we also include convenient plotting utilities for quick visualization:

```python
from scripts.plot_cross_sections import *

# Individual spin group contributions
plot_spin_group_breakdown(compound, 'elastic')

# Compare all reaction types
plot_all_reactions(compound)

# Customizable energy ranges and styling
plot_elastic_comparison(compound, energy_range=(-3, 3), save_path='figure.png')

# But remember - you can always build custom plots using the direct API:
energies = np.logspace(-2, 2, 1000)
elastic_xs = [compound.elasticCrossSection(E) for E in energies]
plt.loglog(energies, elastic_xs, label='Custom Plot')
```

## 🧪 Example Calculations

### Resonance Analysis & Spin Group Physics
```python
# Examine nuclear structure
compound.printSpinGroupInfo()
# → Spin Group 0/2 (0.5, 1): has 3 channels.
# → Spin Group 1/2 (1.5, 1): has 3 channels.

# Calculate at thermal and resonance energies
thermal_energy = 0.0253   # eV (room temperature)
resonance_energy = 0.296  # eV (Pu-239 major resonance)

print("At thermal energy:")
print(f"  Total: {compound.totalCrossSection(thermal_energy):.0f} barns")
print(f"  J=1/2: {compound.spinGroupTotalCrossSection(0, thermal_energy):.0f} barns")
print(f"  J=3/2: {compound.spinGroupTotalCrossSection(1, thermal_energy):.0f} barns")

print(f"\nAt resonance ({resonance_energy} eV):")
print(f"  Total: {compound.totalCrossSection(resonance_energy):.0f} barns")
print(f"  Elastic: {compound.elasticCrossSection(resonance_energy):.0f} barns")
print(f"  Fission: {compound.fissionCrossSection(resonance_energy):.0f} barns")
```

### Energy-dependent Cross Sections & Peak Finding
```python
import numpy as np

# Thermal to fast neutron range
energies = np.logspace(-3, 6, 1000)  # 1 meV to 1 MeV

# Vectorized calculations using our elegant API
total_xs = [compound.totalCrossSection(E) for E in energies]
elastic_xs = [compound.elasticCrossSection(E) for E in energies]
capture_xs = [compound.captureCrossSection(E) for E in energies]

# Physics analysis
peaks = np.where(np.array(total_xs) > 1000)[0]
max_xs = max(total_xs)
max_energy = energies[np.argmax(total_xs)]

print(f"Found {len(peaks)} resonance peaks above 1000 barns")
print(f"Maximum cross section: {max_xs:.0f} barns at {max_energy:.3f} eV")

# Reaction branching ratios
thermal_idx = np.argmin(np.abs(energies - 0.0253))
thermal_total = total_xs[thermal_idx]
print(f"Thermal fission probability: {capture_xs[thermal_idx]/thermal_total:.1%}")
```

### Individual Channel Analysis
```python
# Access spin group details
spin_group_0 = compound.getSpinGroup(0)  # J=1/2 group
print(f"J={spin_group_0.getJ()}, Parity={spin_group_0.getPJ()}")
print(f"Number of channels: {len(spin_group_0.channels())}")
print(f"Number of resonances: {len(spin_group_0.getResonances())}")

# Compare spin group contributions across energy range
energies_thermal = np.logspace(-2, 1, 100)  # 0.01 to 10 eV
j_half_elastic = [compound.spinGroupElasticCrossSection(0, E) for E in energies_thermal]
j_3half_elastic = [compound.spinGroupElasticCrossSection(1, E) for E in energies_thermal]

# Plot spin group interference
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.loglog(energies_thermal, j_half_elastic, label='J=1/2', linewidth=2)
plt.loglog(energies_thermal, j_3half_elastic, label='J=3/2', linewidth=2)
plt.xlabel('Energy (eV)')
plt.ylabel('Elastic Cross Section (barns)')
plt.title('Quantum Mechanical Spin Group Contributions')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
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
