Models For Nuclear Data Evaluation
=======

## Why ?

Models used in nuclear data evaluation, easily callable from Python;

Builds upon solid developments and pedagogical readings :
- SAMMY Manual
- G.Ferran's thesis
- C.Jeannesson's thesis (future devlopment)
- JEFF Report 18 
- ENDF Manual
- Neutron Interaction Theory - A.Foderaro

## Project Structure

- `cpp/`: Contains the C++ code and Pybind11 setup for efficient reconstruction.
- `notebooks/`: Jupyter notebooks for testing and exploring functionalities.

## Setup

### Build Pyrat Extension

Navigate to the main directory and build the extension:

```sh
mkdir build; cd build;
cmake ..
make
```

### Notebook "import pyrat"

```mermaid
  graph TD;
      notebook-->pyrat;
      notebook-->ENDFtk;
      ENDFtk-->Parameters;
      Parameters-->pyrat;
```

## FormalismParametersWrapper Class Diagram

```mermaid
classDiagram
    direction TB
    class FormalismParameters {
        doubles SPI, AP, LAD
        ResonanceParameters parameters
    }

    class MLBWParameters {
        doubles SPI, AP, AWRI
        vector~Llist~ lValues
    }

    class RMatrixParameters {
        vector~JPiList~ JPiGroups
    }

    class Llist {
        double APL, L
        RPlist resonanceParameters
    }

    class RPlist {
        vector~double~ ER, AJ, GN, GG, GFA, GFB
    }

    class JPiList {
        double J, Pi
    }

    FormalismParameters --> ResonanceParameters : parameters
    ResonanceParameters <|-- MLBWParameters
    ResonanceParameters <|-- RMatrixParameters
    MLBWParameters --> Llist : lValues
    Llist --> RPlist : resonanceParameters
    RMatrixParameters --> JPiList : JPiGroups
```


## MLBW Class Diagram


```mermaid
classDiagram
    direction TB
    class MLBW {
        RadialWaveFunctions
        ParticlePair
        MLBWParameters
    }

    class RadialWaveFunctions {
    }

    class ParticlePair {
        double massA, massB
    }

    MLBW --> RadialWaveFunctions : RWF
    MLBW --> ParticlePair : entrancePP
    MLBW --> MLBWParameters : parameters
```

## Notes 

I will stop RML implementation.
>>>>>>> First commit
