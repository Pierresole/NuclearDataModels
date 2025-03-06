// MLBW.cpp
#include "Mlbw.h"

MLBW::MLBW(const MLBWParameters& params) : params_(params), radialWF_(RadialWaveFunctions(params_.AP)), entranceParticlePair_(ParticlePair(params.AWRI)) {

    // Initialize the vectors with the correct size
    P_at_ER.resize(params_.lValues.size());
    S_at_ER.resize(params_.lValues.size());

    // Loop through each L value in MLBWParameters
    for (size_t lIndex = 0; lIndex < params_.lValues.size(); ++lIndex) {
        const auto& lValue = params_.lValues[lIndex];
        double AP = lValue.APL;  // Retrieve the AP value for current L value

        // Initialize vectors for the current l value
        P_at_ER[lIndex].resize(lValue.resonanceParameters.ER.size());
        S_at_ER[lIndex].resize(lValue.resonanceParameters.ER.size());

        // Loop through each resonance within this L value
        for (size_t resIndex = 0; resIndex < lValue.resonanceParameters.ER.size(); ++resIndex) {
            double E_r = lValue.resonanceParameters.ER[resIndex];
            double dk2_Er = 0 ;//entranceParticlePair_.k2(E_r);
            double dRhoTrue2_Er = dk2_Er * std::pow(AP, 2);

            // Compute P and S at E_r
            //radialWF_.P_S_Phi(lValue.L, dRhoTrue2_Er);
            P_at_ER[lIndex][resIndex] = radialWF_.getP(lValue.L);
            S_at_ER[lIndex][resIndex] = radialWF_.getS(lValue.L);
        }
    }
}

double MLBW::radiative_capture_cross_section(double E) {
    double SPI = params_.SPI;
    double dk2 = 0; //entranceParticlePair_.k2(E);
    double sigma_ng = 0.0;
    double dRhoTrue2 = dk2 * std::pow(params_.lValues[0].APL, 2); // Assuming AP from the first L value

    // Precompute P, S, and Phi values once for the current energy E
    for (int l = 0; l <= 4; ++l) {
        radialWF_.P_S_Phi(l, dRhoTrue2);
    }

    // Loop through each L value
    for (size_t lIndex = 0; lIndex < params_.lValues.size(); ++lIndex) {
        const auto& lValue = params_.lValues[lIndex];
        double sigma_l_ng = 0.0;

        // Loop through each resonance within this L value
        for (size_t resIndex = 0; resIndex < lValue.resonanceParameters.ER.size(); ++resIndex) {
            double E_r = lValue.resonanceParameters.ER[resIndex];
            double GN = lValue.resonanceParameters.GN[resIndex];
            double GG = lValue.resonanceParameters.GG[resIndex];

            // Use the precomputed P and S at E_r
            double P_Er = P_at_ER[lIndex][resIndex];
            double S_Er = S_at_ER[lIndex][resIndex];

            // Adjust the neutron width for energy E using the penetration factor
            double widthAtEr = (P_Er / radialWF_.getP(lValue.L)) * GN; // Neutron width scaling

            // Calculate the primed resonance energy using precomputed shift factors
            double dShiftDiff = S_Er - radialWF_.getS(lValue.L);
            double E_r_prime = E_r + (dShiftDiff / (2.0 * P_Er)) * GN;

            double GT_E = widthAtEr + GG;
            double dE = E - E_r_prime;

            // Compute denominator for the cross-section calculation
            double denominator = dE * dE + 0.25 * GT_E * GT_E;

            double GJ = 0.5 * (2 * std::abs(lValue.resonanceParameters.AJ[resIndex]) + 1) / (2 * SPI + 1);

            sigma_l_ng += GJ * (widthAtEr * GG) / denominator;
        }

        // Accumulate contributions for each l state
        sigma_ng += sigma_l_ng;
    }

    // Final cross-section calculation
    return (M_PI / dk2) * sigma_ng;
}

