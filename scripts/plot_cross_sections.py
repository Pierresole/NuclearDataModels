#!/usr/bin/env python3
"""
Cross Section Plotting Utilities for Nuclear Data Models

This module provides convenient plotting functions for visualizing cross sections
computed with the pyRMatrix library. Perfect for quick analysis and demonstrations.

Usage:
    from scripts.plot_cross_sections import *
    
    # Quick plot
    compound = create_compound_from_ReichMoore('file.endf')
    plot_elastic_comparison(compound, save_path='demo.png')
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from pathlib import Path

# Add the build directory to the path for pyRMatrix
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root / 'build' / 'python'))
sys.path.append(str(project_root / 'scripts'))

try:
    import pyRMatrix
    from compoundFromENDFtk import create_compound_from_ReichMoore
except ImportError as e:
    print(f"Warning: Could not import required modules: {e}")
    print("Make sure pyRMatrix is built and compoundFromENDFtk is available")

# Set up matplotlib for publication-quality plots
plt.style.use('seaborn-v0_8' if 'seaborn-v0_8' in plt.style.available else 'default')
COLORS = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#592E83', '#F79D65']

def setup_plot_style():
    """Set up consistent plot styling"""
    plt.rcParams.update({
        'figure.figsize': (10, 6),
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'legend.fontsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'grid.alpha': 0.3,
        'lines.linewidth': 2.5
    })

def plot_elastic_comparison(compound_system, energy_range=(-2, 2), n_points=1000, 
                          save_path=None, show_individual=True, show_total=True):
    """
    Plot elastic cross sections comparing individual spin groups and total.
    
    Args:
        compound_system: CompoundSystem object
        energy_range: tuple of (log10_min, log10_max) in eV
        n_points: number of energy points
        save_path: optional path to save the plot
        show_individual: whether to show individual spin group contributions
        show_total: whether to show total cross section
    
    Returns:
        fig, ax: matplotlib figure and axes objects
    """
    setup_plot_style()
    
    # Create energy array
    energies = np.logspace(energy_range[0], energy_range[1], n_points)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot individual spin groups if requested
    if show_individual:
        try:
            n_spin_groups = len([sg for sg in compound_system.spinGroups()])
            for i in range(min(n_spin_groups, len(COLORS)-1)):  # Reserve one color for total
                try:
                    xs_spingroup = compound_system.spinGroupElasticCrossSection(i, energies.tolist())
                    J_value = compound_system.getSpinGroup(i).getJ()
                    ax.loglog(energies, xs_spingroup, 
                             label=f'Spin Group {i} (J={J_value})', 
                             color=COLORS[i], linewidth=2.5, alpha=0.8)
                except Exception as e:
                    print(f"Warning: Could not plot spin group {i}: {e}")
        except Exception as e:
            print(f"Warning: Could not access spin groups: {e}")
    
    # Plot total if requested
    if show_total:
        try:
            xs_total = compound_system.elasticCrossSection(energies.tolist())
            ax.loglog(energies, xs_total, 
                     label='Total Elastic', linewidth=3, 
                     linestyle='--', color='black')
        except Exception as e:
            print(f"Warning: Could not plot total cross section: {e}")
    
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Elastic Cross Section (barns)')
    ax.set_title('Elastic Cross Section: Individual Spin Groups vs Total')
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    return fig, ax

def plot_all_reactions(compound_system, energy_range=(-2, 2), n_points=1000, save_path=None):
    """
    Plot all reaction types (elastic, capture, fission, total) on the same plot.
    
    Args:
        compound_system: CompoundSystem object
        energy_range: tuple of (log10_min, log10_max) in eV
        n_points: number of energy points
        save_path: optional path to save the plot
    
    Returns:
        fig, ax: matplotlib figure and axes objects
    """
    setup_plot_style()
    
    energies = np.logspace(energy_range[0], energy_range[1], n_points)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    reactions = [
        ('elasticCrossSection', 'Elastic', '#2E86AB'),
        ('captureCrossSection', 'Capture (n,γ)', '#A23B72'),
        ('fissionCrossSection', 'Fission', '#F18F01'),
        ('totalCrossSection', 'Total', '#000000')
    ]
    
    for method_name, label, color in reactions:
        try:
            method = getattr(compound_system, method_name)
            xs = method(energies.tolist())
            linestyle = '--' if label == 'Total' else '-'
            linewidth = 3 if label == 'Total' else 2.5
            ax.loglog(energies, xs, label=label, color=color, 
                     linestyle=linestyle, linewidth=linewidth)
        except Exception as e:
            print(f"Warning: Could not plot {label}: {e}")
    
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Cross Section (barns)')
    ax.set_title('Nuclear Cross Sections by Reaction Type')
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    return fig, ax

def plot_spin_group_breakdown(compound_system, reaction_type='elastic', 
                            energy_range=(-2, 2), n_points=1000, save_path=None):
    """
    Plot individual spin group contributions for a specific reaction type.
    
    Args:
        compound_system: CompoundSystem object
        reaction_type: 'elastic', 'capture', 'fission', or 'total'
        energy_range: tuple of (log10_min, log10_max) in eV
        n_points: number of energy points
        save_path: optional path to save the plot
    
    Returns:
        fig, ax: matplotlib figure and axes objects
    """
    setup_plot_style()
    
    energies = np.logspace(energy_range[0], energy_range[1], n_points)
    method_name = f'spinGroup{reaction_type.capitalize()}CrossSection'
    total_method = f'{reaction_type}CrossSection'
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    try:
        n_spin_groups = len([sg for sg in compound_system.spinGroups()])
        
        # Plot individual spin groups
        for i in range(min(n_spin_groups, len(COLORS)-1)):
            try:
                method = getattr(compound_system, method_name)
                xs_spingroup = method(i, energies.tolist())
                J_value = compound_system.getSpinGroup(i).getJ()
                ax.loglog(energies, xs_spingroup, 
                         label=f'J = {J_value}', 
                         color=COLORS[i], linewidth=2.5)
            except Exception as e:
                print(f"Warning: Could not plot spin group {i}: {e}")
        
        # Plot total
        try:
            total_method_func = getattr(compound_system, total_method)
            xs_total = total_method_func(energies.tolist())
            ax.loglog(energies, xs_total, 
                     label='Total', linewidth=3, 
                     linestyle='--', color='black')
        except Exception as e:
            print(f"Warning: Could not plot total: {e}")
        
    except Exception as e:
        print(f"Error: {e}")
    
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel(f'{reaction_type.capitalize()} Cross Section (barns)')
    ax.set_title(f'{reaction_type.capitalize()} Cross Section by Spin Group')
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    return fig, ax

def quick_demo_plot(endf_file_path, save_path=None):
    """
    Create a quick demonstration plot from an ENDF file.
    Perfect for showcasing the library capabilities.
    
    Args:
        endf_file_path: path to ENDF file
        save_path: optional path to save the plot
    
    Returns:
        fig, ax: matplotlib figure and axes objects
    """
    try:
        # Load compound system
        compound_system = create_compound_from_ReichMoore(endf_file_path)
        
        # Create the demo plot
        fig, ax = plot_elastic_comparison(compound_system, save_path=save_path)
        
        # Add a subtitle with file info
        filename = os.path.basename(endf_file_path)
        ax.text(0.02, 0.98, f'Source: {filename}', 
                transform=ax.transAxes, fontsize=10, 
                verticalalignment='top', alpha=0.7)
        
        return fig, ax
        
    except Exception as e:
        print(f"Error creating demo plot: {e}")
        return None, None

def create_github_demo():
    """
    Create the demonstration plot for GitHub README/research.md
    """
    # Default Pu-239 file path (adjust as needed)
    default_file = '/home/sole-pie01/ndlib/endfb8-neutron/n-094_Pu_239.endf'
    
    if os.path.exists(default_file):
        save_path = project_root / 'images' / 'elastic_cross_section_demo.png'
        save_path.parent.mkdir(exist_ok=True)
        
        print("Creating GitHub demo plot...")
        fig, ax = quick_demo_plot(default_file, str(save_path))
        
        if fig:
            plt.show()
            print(f"Demo plot created and saved to: {save_path}")
        else:
            print("Failed to create demo plot")
    else:
        print(f"Default file not found: {default_file}")
        print("Please provide the path to an ENDF file")

if __name__ == "__main__":
    """
    Command line interface for quick plotting
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot nuclear cross sections')
    parser.add_argument('endf_file', nargs='?', help='Path to ENDF file')
    parser.add_argument('--type', choices=['elastic', 'all', 'breakdown'], 
                       default='elastic', help='Type of plot to create')
    parser.add_argument('--reaction', choices=['elastic', 'capture', 'fission', 'total'], 
                       default='elastic', help='Reaction type for breakdown plot')
    parser.add_argument('--save', help='Path to save the plot')
    parser.add_argument('--demo', action='store_true', 
                       help='Create GitHub demo plot')
    
    args = parser.parse_args()
    
    if args.demo:
        create_github_demo()
    elif args.endf_file:
        try:
            compound_system = create_compound_from_ReichMoore(args.endf_file)
            
            if args.type == 'elastic':
                plot_elastic_comparison(compound_system, save_path=args.save)
            elif args.type == 'all':
                plot_all_reactions(compound_system, save_path=args.save)
            elif args.type == 'breakdown':
                plot_spin_group_breakdown(compound_system, args.reaction, save_path=args.save)
            
            plt.show()
            
        except Exception as e:
            print(f"Error: {e}")
    else:
        print("Please provide an ENDF file or use --demo flag")
        parser.print_help()