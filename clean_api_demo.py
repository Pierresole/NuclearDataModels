#!/usr/bin/env python3
"""
Clean API Showcase Demo - Two focused panels showing the most impressive features
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Add paths
sys.path.append('build/python')
sys.path.append('scripts')

try:
    import pyRMatrix
    from compoundFromENDFtk import create_compound_from_ReichMoore
    print("✅ Successfully imported pyRMatrix modules")
except ImportError as e:
    print(f"❌ Import error: {e}")
    sys.exit(1)

def clean_api_showcase():
    """Create a clean 2-panel showcase of the API"""
    
    # Set up matplotlib for beautiful plots
    plt.style.use('default')
    plt.rcParams.update({
        'figure.figsize': (16, 8),
        'font.size': 13,
        'axes.labelsize': 15,
        'axes.titlesize': 17,
        'legend.fontsize': 13,
        'lines.linewidth': 3,
        'grid.alpha': 0.3
    })
    
    # Colors for different lines
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#592E83', '#4A90A4']
    
    # Try to load nuclear data
    test_files = [
        '/home/sole-pie01/ndlib/endfb8-neutron/n-094_Pu_239.endf',
        '/home/sole-pie01/ndlib/endfb8-neutron/n-029_Cu_063.endf'
    ]
    
    compound_system = None
    filename = "Nuclear Data"
    
    for test_file in test_files:
        if os.path.exists(test_file):
            try:
                print(f"📂 Loading {test_file}...")
                compound_system = create_compound_from_ReichMoore(test_file)
                filename = os.path.basename(test_file).replace('.endf', '').replace('n-', '').replace('_', ' ')
                print(f"✅ Successfully loaded {filename}")
                break
            except Exception as e:
                print(f"⚠️  Failed to load {test_file}: {e}")
                continue
    
    if compound_system is None:
        print("❌ No ENDF files available")
        return
    
    # Create a clean 1x2 subplot layout
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Energy range
    energies = np.logspace(-2, 2, 1000)     # 0.01 to 100 eV
    
    # Left Panel: Multiple Reaction Types - Shows API elegance
    print("🧮 Computing cross sections by reaction type...")
    try:
        elastic_xs = [compound_system.elasticCrossSection(E) for E in energies]
        capture_xs = [compound_system.captureCrossSection(E) for E in energies]
        fission_xs = [compound_system.fissionCrossSection(E) for E in energies]
        total_xs = [compound_system.totalCrossSection(E) for E in energies]
        
        ax1.loglog(energies, elastic_xs, label='Elastic Scattering', color=colors[0], linewidth=3)
        ax1.loglog(energies, capture_xs, label='Radiative Capture (n,γ)', color=colors[1], linewidth=3)
        ax1.loglog(energies, fission_xs, label='Fission', color=colors[2], linewidth=3)
        ax1.loglog(energies, total_xs, label='Total', color='black', linewidth=4, linestyle='--', alpha=0.8)
        
        ax1.set_xlabel('Energy (eV)', fontsize=15)
        ax1.set_ylabel('Cross Section (barns)', fontsize=15)
        ax1.set_title('Nuclear Reaction Types\ncompound.elasticCrossSection(E)', fontsize=17, pad=20)
        ax1.legend(fontsize=13, loc='upper right')
        ax1.grid(True, alpha=0.3)
        
        # Add API showcase text
        api_text = ("# Elegant API calls:\n"
                   "elastic = compound.elasticCrossSection(E)\n"
                   "capture = compound.captureCrossSection(E)\n"
                   "fission = compound.fissionCrossSection(E)\n"
                   "total = compound.totalCrossSection(E)")
        
        ax1.text(0.02, 0.02, api_text, transform=ax1.transAxes, fontsize=11,
                verticalalignment='bottom', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.9, edgecolor='gray'))
        
        print(f"   ✓ Computed {len(energies)} energy points across 4 reaction types")
        
    except Exception as e:
        print(f"   ⚠️  Error computing cross sections: {e}")
    
    # Right Panel: Spin Group Decomposition - Shows quantum physics
    print("🌟 Computing spin group contributions...")
    try:
        energies_thermal = np.logspace(-2, 1, 500)   # Focus on thermal/epithermal region
        n_spin_groups = len([sg for sg in compound_system.spinGroups()])
        
        for i in range(min(n_spin_groups, 4)):  # Limit to first 4 spin groups
            try:
                sg = compound_system.getSpinGroup(i)
                J_value = sg.getJ()
                parity = sg.getPJ()
                
                sg_elastic = [compound_system.spinGroupElasticCrossSection(i, E) for E in energies_thermal]
                
                parity_symbol = "+" if parity > 0 else "-"
                ax2.loglog(energies_thermal, sg_elastic, 
                          label=f'J = {J_value}{parity_symbol}', 
                          color=colors[i], linewidth=3)
                
                print(f"   ✓ Spin group {i}: J={J_value}, Parity={parity}")
                
            except Exception as e:
                print(f"   ⚠️  Error with spin group {i}: {e}")
        
        # Total elastic for comparison
        total_elastic = [compound_system.elasticCrossSection(E) for E in energies_thermal]
        ax2.loglog(energies_thermal, total_elastic, 
                  label='Total Elastic', color='black', 
                  linewidth=4, linestyle='--', alpha=0.8)
        
        ax2.set_xlabel('Energy (eV)', fontsize=15)
        ax2.set_ylabel('Elastic Cross Section (barns)', fontsize=15)
        ax2.set_title('Spin Group Decomposition (Quantum Physics)\ncompound.spinGroupElasticCrossSection(i, E)', fontsize=17, pad=20)
        ax2.legend(fontsize=13)
        ax2.grid(True, alpha=0.3)
        
        # Add quantum physics showcase text
        quantum_text = ("# Quantum mechanics access:\n"
                       "for i, sg in enumerate(compound.spinGroups()):\n"
                       "    J = compound.getSpinGroup(i).getJ()\n"
                       "    xs = compound.spinGroupElasticCrossSection(i, E)")
        
        ax2.text(0.02, 0.02, quantum_text, transform=ax2.transAxes, fontsize=11,
                verticalalignment='bottom', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.9, edgecolor='gray'))
        
    except Exception as e:
        print(f"   ⚠️  Error computing spin groups: {e}")
    
    # Add main title
    fig.suptitle(f'pyRMatrix API Showcase: {filename.title()} Nuclear Data', fontsize=20, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Save the plot
    save_path = 'images/api_showcase_demo.png'
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"🎨 Clean API showcase plot saved to: {save_path}")
    
    # Also save as SVG
    save_path_svg = 'images/api_showcase_demo.svg'
    plt.savefig(save_path_svg, bbox_inches='tight', facecolor='white')
    print(f"🎨 Clean API showcase plot saved to: {save_path_svg}")
    
    plt.show()

if __name__ == "__main__":
    print("🚀 Creating Clean pyRMatrix API Showcase")
    print("=" * 50)
    clean_api_showcase()
    print("=" * 50)
    print("✅ Clean API Showcase Complete!")
