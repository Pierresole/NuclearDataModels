#!/usr/bin/env python3
"""
Create a demonstration plot for the README
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from pathlib import Path

# Add paths
sys.path.append('build/python')
sys.path.append('scripts')

try:
    import pyRMatrix
    from compoundFromENDFtk import create_compound_from_ReichMoore
    print("Successfully imported modules")
except ImportError as e:
    print(f"Import error: {e}")
    sys.exit(1)

# Set up matplotlib for nice plots
plt.style.use('default')
plt.rcParams.update({
    'figure.figsize': (12, 8),
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'legend.fontsize': 12,
    'lines.linewidth': 2.5,
    'grid.alpha': 0.3
})

def create_demo_plot():
    """Create demonstration plot"""
    
    # Colors for different lines
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#592E83']
    
    # Try a simple test file or create synthetic data if ENDF file isn't available
    test_files = [
        '/home/sole-pie01/ndlib/endfb8-neutron/n-094_Pu_239.endf',
        '/home/sole-pie01/ndlib/endfb8-neutron/n-029_Cu_063.endf'
    ]
    
    compound_system = None
    filename = "Synthetic Data"
    
    for test_file in test_files:
        if os.path.exists(test_file):
            try:
                print(f"Loading {test_file}...")
                compound_system = create_compound_from_ReichMoore(test_file)
                filename = os.path.basename(test_file)
                print(f"Successfully loaded {filename}")
                break
            except Exception as e:
                print(f"Failed to load {test_file}: {e}")
                continue
    
    if compound_system is None:
        print("No ENDF files available, creating synthetic demonstration plot...")
        create_synthetic_plot()
        return
    
    # Create energy array
    energies = np.logspace(-2, 2, 1000)  # 0.01 eV to 100 eV
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    try:
        # Plot different cross section types
        reactions = [
            ('elasticCrossSection', 'Elastic Scattering', colors[0]),
            ('captureCrossSection', 'Radiative Capture (n,γ)', colors[1]),
            ('fissionCrossSection', 'Fission', colors[2]),
            ('totalCrossSection', 'Total', 'black')
        ]
        
        for method_name, label, color in reactions:
            try:
                method = getattr(compound_system, method_name)
                xs_values = []
                for E in energies:
                    try:
                        xs_values.append(method(E))
                    except:
                        xs_values.append(0.0)
                
                # Only plot if we have meaningful data
                if max(xs_values) > 1e-10:
                    linestyle = '--' if label == 'Total' else '-'
                    linewidth = 3 if label == 'Total' else 2.5
                    ax.loglog(energies, xs_values, label=label, color=color, 
                             linestyle=linestyle, linewidth=linewidth)
                    print(f"Plotted {label}: max = {max(xs_values):.2e} barns")
                
            except Exception as e:
                print(f"Warning: Could not plot {label}: {e}")
    
        ax.set_xlabel('Energy (eV)', fontsize=14)
        ax.set_ylabel('Cross Section (barns)', fontsize=14)
        ax.set_title(f'Nuclear Cross Sections - {filename}', fontsize=16)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=12)
        
        # Add annotations
        ax.text(0.02, 0.98, 'pyRMatrix Nuclear Data Library', 
                transform=ax.transAxes, fontsize=10, 
                verticalalignment='top', alpha=0.7,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        
        plt.tight_layout()
        
        # Save the plot
        save_path = 'images/cross_section_demo.png'
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"Demo plot saved to: {save_path}")
        
        # Also save as SVG for better quality
        save_path_svg = 'images/cross_section_demo.svg'
        plt.savefig(save_path_svg, bbox_inches='tight', facecolor='white')
        print(f"Demo plot saved to: {save_path_svg}")
        
        plt.show()
        
    except Exception as e:
        print(f"Error creating plot: {e}")
        create_synthetic_plot()

def create_synthetic_plot():
    """Create a synthetic demonstration plot if no ENDF data is available"""
    print("Creating synthetic demonstration plot...")
    
    # Create synthetic cross section data that looks realistic
    energies = np.logspace(-2, 2, 1000)
    
    # Synthetic elastic cross section with 1/v behavior and resonances  
    elastic = 10 / np.sqrt(energies) + 5 * np.exp(-((np.log10(energies) - 0.5) / 0.3)**2) + 50 * np.exp(-((np.log10(energies) - 1.2) / 0.2)**2)
    
    # Synthetic capture with 1/v and resonances
    capture = 50 / np.sqrt(energies) + 200 * np.exp(-((np.log10(energies) - 0.3) / 0.2)**2)
    
    # Synthetic fission (only for fissile materials, above threshold)
    fission = np.where(energies > 0.1, 5 / np.sqrt(energies) + 100 * np.exp(-((np.log10(energies) - 1.0) / 0.4)**2), 0)
    
    # Total is sum
    total = elastic + capture + fission
    
    colors = ['#2E86AB', '#A23B72', '#F18F01', 'black']
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.loglog(energies, elastic, label='Elastic Scattering', color=colors[0], linewidth=2.5)
    ax.loglog(energies, capture, label='Radiative Capture (n,γ)', color=colors[1], linewidth=2.5)
    ax.loglog(energies, fission, label='Fission', color=colors[2], linewidth=2.5)
    ax.loglog(energies, total, label='Total', color=colors[3], linewidth=3, linestyle='--')
    
    ax.set_xlabel('Energy (eV)', fontsize=14)
    ax.set_ylabel('Cross Section (barns)', fontsize=14)
    ax.set_title('Nuclear Cross Sections (Demonstration)', fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)
    
    # Add annotations
    ax.text(0.02, 0.98, 'pyRMatrix Nuclear Data Library\n(Synthetic Demonstration Data)', 
            transform=ax.transAxes, fontsize=10, 
            verticalalignment='top', alpha=0.7,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    plt.tight_layout()
    
    # Save the plot
    save_path = 'images/cross_section_demo.png'
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Demo plot saved to: {save_path}")
    
    # Also save as SVG
    save_path_svg = 'images/cross_section_demo.svg'
    plt.savefig(save_path_svg, bbox_inches='tight', facecolor='white')
    print(f"Demo plot saved to: {save_path_svg}")
    
    plt.show()

if __name__ == "__main__":
    create_demo_plot()
