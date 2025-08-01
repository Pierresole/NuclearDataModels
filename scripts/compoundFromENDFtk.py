try:
    import ENDFtk
except ImportError:
    raise ImportError("The ENDFtk module is required but not installed.")
import numpy as np
import sys
sys.path.append('/home/sole-pie01/codes/NuclearDataModels/build/python')
import pyRMatrix

def create_compound_from_RMatrix(tape_string):
    # Parse the ENDF tape from the input string
    tape = ENDFtk.tree.Tape.from_file(tape_string)
    mat_number = tape.material_numbers[0]
    mf2 = tape.MAT(mat_number).MF(2).MT(151).parse()
    resonance_range = mf2.isotopes[0].resonance_ranges[0]
    
        
    # Find the index of 2 in the list of MT values
    indexEntrancePP = resonance_range.parameters.particle_pairs.MT[:].index(2)
    IA = resonance_range.parameters.particle_pairs.IA[indexEntrancePP]
    IB = resonance_range.parameters.particle_pairs.IB[indexEntrancePP]

    # Create entrance particle pair (neutron)
    entrance_pp = pyRMatrix.ParticlePair(
        resonance_range.parameters.particle_pairs.MA[indexEntrancePP],  # mass1
        resonance_range.parameters.particle_pairs.MB[indexEntrancePP],  # mass2
        np.abs(resonance_range.parameters.particle_pairs.IA[indexEntrancePP]),  # spin1
        np.abs(resonance_range.parameters.particle_pairs.IB[indexEntrancePP]),  # spin2
        resonance_range.parameters.particle_pairs.Q[indexEntrancePP],   # Q value
        int(resonance_range.parameters.particle_pairs.PA[indexEntrancePP]) if IA == 0 else int(np.sign(IA)),  # parity1
        int(resonance_range.parameters.particle_pairs.PB[indexEntrancePP]) if IB == 0 else int(np.sign(IB)),  # parity2
        resonance_range.parameters.particle_pairs.MT[indexEntrancePP]   # MT
    )


    
    # Create SpinGroups
    spin_groups = []
    for iJP in range(resonance_range.parameters.NJS):
        AJ = float(resonance_range.parameters.spin_groups[iJP].AJ)
        PJ = int(resonance_range.parameters.spin_groups[iJP].PJ) if AJ == 0 else np.sign(AJ)
        sg = pyRMatrix.SpinGroup(np.abs(AJ), int(PJ))
        
        # Add channels
        for iCH in range(resonance_range.parameters.spin_groups[iJP].NCH):
            ppi_idx = resonance_range.parameters.spin_groups[iJP].channels.PPI[iCH] - 1
            pp = pyRMatrix.ParticlePair(
                float(resonance_range.parameters.particle_pairs.MA[ppi_idx]),
                float(resonance_range.parameters.particle_pairs.MB[ppi_idx]),
                float(resonance_range.parameters.particle_pairs.IA[ppi_idx]),
                float(resonance_range.parameters.particle_pairs.IB[ppi_idx]),
                float(resonance_range.parameters.particle_pairs.Q[ppi_idx]),
                int(resonance_range.parameters.particle_pairs.PA[ppi_idx]),
                int(resonance_range.parameters.particle_pairs.PB[ppi_idx]),
                int(resonance_range.parameters.particle_pairs.MT[ppi_idx])
            )
            
            ch = pyRMatrix.Channel(
                pp,
                resonance_range.parameters.spin_groups[iJP].channels.L[iCH],
                resonance_range.parameters.spin_groups[iJP].channels.APE[iCH],
                resonance_range.parameters.spin_groups[iJP].channels.APT[iCH],
                resonance_range.parameters.spin_groups[iJP].channels.SCH[iCH]
            )
            isEliminated = int(resonance_range.parameters.particle_pairs.MT[ppi_idx]) == 102
            sg.addChannel(ch, isEliminated)
        
        # Add resonances
        for iREZ in range(resonance_range.parameters.spin_groups[iJP].NRS):
            Er = resonance_range.parameters.spin_groups[iJP].parameters.ER[iREZ]
            vGAM = [g for g in resonance_range.parameters.spin_groups[iJP].parameters.GAM[iREZ][:resonance_range.parameters.spin_groups[iJP].NCH]]
            sg.addResonance(Er, vGAM, resonance_range.parameters.IFG)
        
        spin_groups.append(sg)
    
    # Create the CompoundNucleus object
    compound_system = pyRMatrix.CompoundSystem(entrance_pp, spin_groups)
    
    return compound_system


def create_compound_from_ReichMoore(tape_string):
    # Parse the ENDF tape from the input string
    tape = ENDFtk.tree.Tape.from_file(tape_string)
    mat_number = tape.material_numbers[0]
    mf2 = tape.MAT(mat_number).MF(2).MT(151).parse()
    resonance_range = mf2.isotopes[0].resonance_ranges[0]
    
    dMassTarget = mf2.AWR
    dSpinTarget = resonance_range.parameters.SPI
    iParityTarget = 1
    iParityNeutron = 1
    # Create entrance particle pair (neutron)
    entrance_pp = pyRMatrix.ParticlePair.neutron_incident(
        mass2 = dMassTarget,  # mass2
        spin2 = dSpinTarget,  # spin2
        Q = 0,               # because entrance channel
        parity2 = iParityTarget,        
        MT = 2               # because entrance channel
    )

    capture_pp = pyRMatrix.ParticlePair(
        mass1 = 0.0,  # mass1
        mass2 = dMassTarget,  # mass2
        spin1 = 0.0,  # spin1
        spin2 = dSpinTarget,  # spin2
        Q = 0.0,      # Q value
        parity1 = 1,  # parity of the neutron
        parity2 = iParityTarget,  
        MT = 102      # fission channel
    )

    fission1_pp = pyRMatrix.ParticlePair(
        mass1 = dMassTarget/2,
        mass2 = dMassTarget/2,
        spin1 = dSpinTarget,
        spin2 = 0,
        Q = 0.0,      # Q value
        parity1 = iParityTarget,  # parity of the neutron
        parity2 = 1,  # does not matter
        MT = 18      # fission channel
    )

    fission2_pp = pyRMatrix.ParticlePair(
        mass1 = dMassTarget/2,
        mass2 = dMassTarget/2,
        spin1 = dSpinTarget,
        spin2 = 0,
        Q = 0.0,      # Q value
        parity1 = iParityTarget,
        parity2 = 1,  # does not matter
        MT = 18      # fission channel
    )

    #  Returns the channel spin s associated with a value AJ, an orbital angular
    #       momentum l and a target nucleus spin I
    # 
    # It is useful for the formalisms LRU=1 and LRF=1, 2 or 3 of ENDF-6 format
    # In those formalism, AJ is the J-value with a sign (used to determine the spin)
    # 
    # The basis of this method is that in these formalisms, channel spin is supposed
    #   to be conserved, and thus we can consider that the spin we are looking for
    #   is the spin of an entrance channel neutron (spin 1/2) + target nucleus
    # 
    # Therefore the spin is s = |I-0.5| or I+0.5
    # 
    # And we have the vectorial sum : J = s + l
    def CalculatesSpin(J, L, I):
        if(I==0):
            return 0.5

        if(L==0): 
            return np.abs(J)

        if(J==0):
            return L

        if(np.abs(np.abs(J)-L) >= I): 
            return I+0.5

        if(np.abs(J)+L <= I):
            return np.abs(I-0.5)

        if(J<0):
            return np.abs(I-0.5)
        else:
            return I+0.5
    
    from collections import defaultdict

    # Step 1: Group resonances by (AJ, PJ)
    spin_group_dict = defaultdict(list)

    for ilg, lvalue in enumerate(resonance_range.parameters.l_values.to_list()):
        for iResonance in range(lvalue.NRS):
            ER = lvalue.ER[iResonance]
            AJ = lvalue.spin_values[iResonance]
            GN = lvalue.neutron_widths[iResonance]
            GG = lvalue.gamma_widths[iResonance]
            GFA = lvalue.first_fission_widths[iResonance]
            GFB = lvalue.second_fission_widths[iResonance]
            PJ = ((-1)**lvalue.L) * iParityTarget * iParityNeutron
            ChannelSpin = CalculatesSpin(AJ, lvalue.L, dSpinTarget)
            # Store all info needed for resonance
            spin_group_dict[(AJ, PJ)].append({
                'ER': ER,
                'L': lvalue.L,
                'APL': lvalue.APL,
                'ChannelSpin': ChannelSpin,
                'widths': [GG, GN, GFA, GFB]
            })

    # Step 2: Build spin groups
    spin_groups = []
    for (AJ, PJ), resonances in spin_group_dict.items():
        sg = pyRMatrix.SpinGroup(np.abs(AJ), int(PJ))
        # Add 4 channels (using the first resonance's L/APL/ChannelSpin as example)
        chs = []
        for pp, pp_obj in zip(
            [capture_pp, entrance_pp, fission1_pp, fission2_pp],
            ['capture', 'entrance', 'fission1', 'fission2']
        ):
            ch = pyRMatrix.Channel(
                particle_pair=pp,
                l=resonances[0]['L'],
                effective_radius=resonances[0]['APL'],
                true_radius=resonances[0]['APL'],
                channel_spin=resonances[0]['ChannelSpin']
            )
            # For capture channel, isEliminated is True if MT==102
            isEliminated = (pp.MT() == 102)
            sg.addChannel(ch, isEliminated)
            chs.append(ch)
        # Add all resonances to this spin group
        for res in resonances:
            sg.addResonance(res['ER'], res['widths'], isReduced=False)
        spin_groups.append(sg)

    compound_system = pyRMatrix.CompoundSystem(entrance_pp, spin_groups)
    return compound_system