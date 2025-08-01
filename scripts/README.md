# Nuclear Data Models - Scripts

This directory contains utility scripts for working with the Nuclear Data Models library.

## plot_cross_sections.py

A comprehensive plotting utility for visualizing nuclear cross sections computed with pyRMatrix.

### Quick Start

```python
from scripts.plot_cross_sections import *

# Load your data
compound_system = create_compound_from_ReichMoore('your_file.endf')

# Create plots with one line
plot_elastic_comparison(compound_system)          # Individual + total elastic
plot_all_reactions(compound_system)               # All reaction types
plot_spin_group_breakdown(compound_system, 'fission')  # Fission by spin group
```

### Command Line Usage

```bash
# Create demo plot for GitHub
python scripts/plot_cross_sections.py --demo

# Plot from ENDF file
python scripts/plot_cross_sections.py your_file.endf --type elastic --save demo.png
python scripts/plot_cross_sections.py your_file.endf --type all
python scripts/plot_cross_sections.py your_file.endf --type breakdown --reaction fission
```

### Functions

- **`plot_elastic_comparison()`**: Compare individual spin groups vs total elastic cross section
- **`plot_all_reactions()`**: Show elastic, capture, fission, and total cross sections
- **`plot_spin_group_breakdown()`**: Individual spin group contributions for any reaction type
- **`quick_demo_plot()`**: Generate a quick demo from an ENDF file
- **`create_github_demo()`**: Create the standard demo plot for documentation

### Features

- 📊 Publication-quality plots with consistent styling
- 🎨 Automatic color schemes and professional formatting
- 💾 High-resolution PNG export for papers/presentations
- 🔧 Customizable energy ranges and plot parameters
- ⚡ Error handling for missing data or invalid files

### GitHub Demo

To create the demo plot used in the research.md:

```python
from scripts.plot_cross_sections import create_github_demo
create_github_demo()
```

This generates `images/elastic_cross_section_demo.png` perfect for your GitHub page.

## compoundFromENDFtk.py

Utility functions for loading nuclear data from ENDF files into pyRMatrix compound systems.

### Usage

```python
from scripts.compoundFromENDFtk import create_compound_from_ReichMoore

# Load Reich-Moore resonance parameters
compound = create_compound_from_ReichMoore('/path/to/endf/file.endf')

# Use with pyRMatrix
cross_sections = compound.elasticCrossSection([0.1, 1.0, 10.0])
```