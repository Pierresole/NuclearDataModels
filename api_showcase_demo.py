#!/usr/bin/env python3
"""
API Showcase Demo - Demonstrates the elegant nuclear physics API
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

def api_showcase_demo():
    """Demonstrate the elegant compound system API"""
    
    # Set up matplotlib for beautiful plots
    plt.style.use('default')
    plt.rcParams.update({
        'figure.figsize': (14, 10),
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'legend.fontsize': 12,
        'lines.linewidth': 2.5,
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
    
    # Create a 2x2 subplot layout showing API elegance
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Energy ranges
    energies_thermal = np.logspace(-2, 1, 500)   # 0.01 to 10 eV (thermal/epithermal)
    energies_wide = np.logspace(-2, 2, 1000)     # 0.01 to 100 eV
    
    # 1. Cross Section Types (top-left)
    print("🧮 Computing cross sections by reaction type...")
    try:
        elastic_xs = [compound_system.elasticCrossSection(E) for E in energies_wide]
        capture_xs = [compound_system.captureCrossSection(E) for E in energies_wide]
        fission_xs = [compound_system.fissionCrossSection(E) for E in energies_wide]
        total_xs = [compound_system.totalCrossSection(E) for E in energies_wide]
        
        ax1.loglog(energies_wide, elastic_xs, label='Elastic Scattering', color=colors[0], linewidth=2.5)
        ax1.loglog(energies_wide, capture_xs, label='Radiative Capture (n,γ)', color=colors[1], linewidth=2.5)
        ax1.loglog(energies_wide, fission_xs, label='Fission', color=colors[2], linewidth=2.5)
        ax1.loglog(energies_wide, total_xs, label='Total', color='black', linewidth=3, linestyle='--')
        
        ax1.set_xlabel('Energy (eV)')
        ax1.set_ylabel('Cross Section (barns)')
        ax1.set_title('API Demo: Multiple Reaction Types\ncompound.elasticCrossSection(E)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        print(f"   ✓ Computed {len(energies_wide)} energy points across 4 reaction types")
        
    except Exception as e:
        print(f"   ⚠️  Error computing cross sections: {e}")
    
    # 2. Spin Group Decomposition (top-right)
    print("🌟 Computing spin group contributions...")
    try:
        n_spin_groups = len([sg for sg in compound_system.spinGroups()])
        spin_group_elastic = []
        
        for i in range(min(n_spin_groups, 4)):  # Limit to first 4 spin groups
            try:
                sg = compound_system.getSpinGroup(i)
                J_value = sg.getJ()
                parity = sg.getPJ()
                
                sg_elastic = [compound_system.spinGroupElasticCrossSection(i, E) for E in energies_thermal]
                
                parity_symbol = "+" if parity > 0 else "-"
                ax2.loglog(energies_thermal, sg_elastic, 
                          label=f'J={J_value}{parity_symbol}', 
                          color=colors[i], linewidth=2.5)
                spin_group_elastic.append(sg_elastic)
                
                print(f"   ✓ Spin group {i}: J={J_value}, Parity={parity}")
                
            except Exception as e:
                print(f"   ⚠️  Error with spin group {i}: {e}")
        
        # Total elastic for comparison
        total_elastic = [compound_system.elasticCrossSection(E) for E in energies_thermal]
        ax2.loglog(energies_thermal, total_elastic, 
                  label='Total Elastic', color='black', 
                  linewidth=3, linestyle='--', alpha=0.8)
        
        ax2.set_xlabel('Energy (eV)')
        ax2.set_ylabel('Elastic Cross Section (barns)')
        ax2.set_title('API Demo: Spin Group Decomposition\ncompound.spinGroupElasticCrossSection(i, E)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
    except Exception as e:
        print(f"   ⚠️  Error computing spin groups: {e}")
    
    # 3. Energy Point Analysis (bottom-left)
    print("📊 Computing point-wise comparisons...")
    try:
        # Select specific energies for demonstration
        demo_energies = [0.0253, 0.1, 1.0, 10.0]  # thermal, epithermal, intermediate, fast
        energy_labels = ['Thermal\n(0.025 eV)', 'Epithermal\n(0.1 eV)', 'Intermediate\n(1 eV)', 'Fast\n(10 eV)']
        
        # Get cross sections at these energies
        elastic_points = [compound_system.elasticCrossSection(E) for E in demo_energies]
        capture_points = [compound_system.captureCrossSection(E) for E in demo_energies]
        fission_points = [compound_system.fissionCrossSection(E) for E in demo_energies]
        
        x = np.arange(len(demo_energies))
        width = 0.25
        
        ax3.bar(x - width, elastic_points, width, label='Elastic', color=colors[0], alpha=0.8)
        ax3.bar(x, capture_points, width, label='Capture', color=colors[1], alpha=0.8)
        ax3.bar(x + width, fission_points, width, label='Fission', color=colors[2], alpha=0.8)
        
        ax3.set_xlabel('Neutron Energy Regime')
        ax3.set_ylabel('Cross Section (barns)')
        ax3.set_title('API Demo: Point-wise Calculations\ncompound.elasticCrossSection(energy)')
        ax3.set_xticks(x)
        ax3.set_xticklabels(energy_labels)
        ax3.set_yscale('log')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        print(f"   ✓ Computed cross sections at {len(demo_energies)} characteristic energies")
        
    except Exception as e:
        print(f"   ⚠️  Error computing energy points: {e}")
    
    # 4. Nuclear Structure Info (bottom-right) - Text summary
    print("📋 Extracting nuclear structure information...")
    try:
        ax4.axis('off')  # Turn off axes for text display
        
        # Collect system information
        info_text = f"Nuclear System: {filename}\n\n"
        
        # Entrance particle pair info
        entrance_pp = compound_system.entranceParticlePair()
        info_text += f"Entrance Channel:\n"
        info_text += f"  Projectile mass: {entrance_pp.mass1():.3f} amu\n"
        info_text += f"  Target mass: {entrance_pp.mass2():.1f} amu\n"
        info_text += f"  Projectile spin: {entrance_pp.spin1()}\n"
        info_text += f"  Target spin: {entrance_pp.spin2()}\n\n"
        
        # Spin group summary
        info_text += f"Spin Groups:\n"
        n_groups = len([sg for sg in compound_system.spinGroups()])
        total_channels = 0
        total_resonances = 0
        
        for i in range(min(n_groups, 5)):  # Show first 5 groups
            try:
                sg = compound_system.getSpinGroup(i)
                J = sg.getJ()
                parity = sg.getPJ()
                n_channels = len(sg.channels())
                n_resonances = len(sg.getResonances())
                
                parity_symbol = "+" if parity > 0 else "-"
                info_text += f"  J={J}{parity_symbol}: {n_channels} channels, {n_resonances} resonances\n"
                
                total_channels += n_channels
                total_resonances += n_resonances
                
            except Exception as e:
                info_text += f"  Group {i}: Error accessing data\n"
        
        if n_groups > 5:
            info_text += f"  ... and {n_groups - 5} more groups\n"
        
        info_text += f"\nTotal: {total_channels} channels, {total_resonances} resonances\n\n"
        
        # Sample calculations
        thermal_E = 0.0253
        info_text += f"Sample Calculations (E = {thermal_E} eV):\n"
        try:
            elastic = compound_system.elasticCrossSection(thermal_E)
            capture = compound_system.captureCrossSection(thermal_E)
            fission = compound_system.fissionCrossSection(thermal_E)
            total = compound_system.totalCrossSection(thermal_E)
            
            info_text += f"  σ_elastic = {elastic:.1f} barns\n"
            info_text += f"  σ_capture = {capture:.1f} barns\n"
            info_text += f"  σ_fission = {fission:.1f} barns\n"
            info_text += f"  σ_total = {total:.1f} barns\n\n"
            
        except Exception as e:
            info_text += f"  Error in calculations: {e}\n\n"
        
        info_text += "API Methods Demonstrated:\n"
        info_text += "• compound.elasticCrossSection(E)\n"
        info_text += "• compound.captureCrossSection(E)\n"
        info_text += "• compound.fissionCrossSection(E)\n"
        info_text += "• compound.spinGroupElasticCrossSection(i, E)\n"
        info_text += "• compound.getSpinGroup(i).getJ()\n"
        info_text += "• compound.entranceParticlePair()\n"
        info_text += "• compound.spinGroups()\n"
        
        ax4.text(0.05, 0.95, info_text, transform=ax4.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
        
        ax4.set_title('API Demo: Nuclear Structure Access\ncompound.getSpinGroup(i).getJ()', fontsize=14)
        
        print(f"   ✓ Extracted information from {n_groups} spin groups")
        
    except Exception as e:
        print(f"   ⚠️  Error extracting nuclear structure: {e}")
    
    # Add main title and save
    fig.suptitle('pyRMatrix API Showcase: Elegant Nuclear Physics Calculations', fontsize=18, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    
    # Save the plot
    save_path = 'images/api_showcase_demo.png'
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"🎨 API showcase plot saved to: {save_path}")
    
    # Also save as SVG
    save_path_svg = 'images/api_showcase_demo.svg'
    plt.savefig(save_path_svg, bbox_inches='tight', facecolor='white')
    print(f"🎨 API showcase plot saved to: {save_path_svg}")
    
    plt.show()

if __name__ == "__main__":
    print("🚀 Starting pyRMatrix API Showcase Demo")
    print("=" * 50)
    api_showcase_demo()
    print("=" * 50)
    print("✅ API Showcase Demo Complete!")
