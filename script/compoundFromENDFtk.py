try:
    import ENDFtk
except ImportError:
    raise ImportError("The ENDFtk module is required but not installed.")

import sys
sys.path.append('/home/sole-pie01/mycodes/pyrat/build/python')
import pyrat

def create_compound_nucleus(tape_string):
    # Parse the ENDF tape from the input string
    tape = ENDFtk.tree.Tape.from_file(tape_string)
    mat_number = tape.material_numbers[0]
    mf2 = tape.MAT(mat_number).MF(2).MT(151).parse()
    resonance_range = mf2.isotopes[0].resonance_ranges[0]
    
    # Prepare ParticlePair instances
    particle_pairs = []
    for ipp in range(resonance_range.parameters.particle_pairs.NPP):
        particle_pair = pyrat.ParticlePair(
            float(resonance_range.parameters.particle_pairs.MA[ipp]),
            float(resonance_range.parameters.particle_pairs.MB[ipp]),
            float(resonance_range.parameters.particle_pairs.IA[ipp]),
            float(resonance_range.parameters.particle_pairs.IB[ipp]),
            float(resonance_range.parameters.particle_pairs.Q[ipp]),
            int(resonance_range.parameters.particle_pairs.PA[ipp]),
            int(resonance_range.parameters.particle_pairs.PB[ipp]),
            int(resonance_range.parameters.particle_pairs.MT[ipp])
        )
        particle_pairs.append(particle_pair)
    
    # Prepare SpinGroup instances
    spin_groups = []
    for iJP in range(resonance_range.parameters.NJS):

        # Collect channel data
        my_channels = []
        for iCH in range(resonance_range.parameters.spin_groups[iJP].NCH):
            channel = pyrat.Channel(
                particle_pairs[resonance_range.parameters.spin_groups[iJP].channels.PPI[iCH] - 1],
                resonance_range.parameters.spin_groups[iJP].channels.L[iCH],
                resonance_range.parameters.spin_groups[iJP].channels.APE[iCH],
                resonance_range.parameters.spin_groups[iJP].channels.APT[iCH],
                resonance_range.parameters.spin_groups[iJP].channels.SCH[iCH]
            )
            my_channels.append(channel)
            # spin_group.addChannel(channel)

        # Collect resonance data
        my_resonances = []
        for iREZ in range(resonance_range.parameters.spin_groups[iJP].NRS):
            resonance = pyrat.Resonance(
                resonance_range.parameters.spin_groups[iJP].parameters.ER[iREZ],
                resonance_range.parameters.spin_groups[iJP].parameters.GAM[iREZ]
            )
            my_resonances.append(resonance)
            # spin_group.addResonance(resonance)

        spin_group = pyrat.SpinGroup(
            J=float(resonance_range.parameters.spin_groups[iJP].AJ),
            PJ=int(resonance_range.parameters.spin_groups[iJP].PJ),
            channels=my_channels,
            resonances=my_resonances
        )
        
        spin_groups.append(spin_group)
    
    # Create the CompoundNucleus object
    entrance_particle_pair = particle_pairs[1]  # Assuming the first particle pair is the entrance channel
    print(entrance_particle_pair.MT())
    compound_nucleus = pyrat.CompoundSystem(entrance_particle_pair, spin_groups)
    # for spin_group in spin_groups:
    #     compound_nucleus.addSpinGroup(spin_group)
    
    return compound_nucleus