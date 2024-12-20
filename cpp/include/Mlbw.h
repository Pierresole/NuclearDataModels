// MLBW.h
#ifndef MLBW_H
#define MLBW_H

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "FormalismParameters.h"   
#include "RadialWaveFunctions.h"
#include "ParticlePair.h"

class MLBW {
public:
    MLBW(const MLBWParameters& params);

    // Methods to compute cross-sections
    double radiative_capture_cross_section(double E);

private:

    double APTrue_; // Scattering radius
    double APEff_;

    // Precomputed values
    std::vector<std::vector<double>> P_at_ER; // Penetration factors at E_r for l and i-th reso
    std::vector<std::vector<double>> S_at_ER; // Shift factors at E_r for l and i-th reso

    RadialWaveFunctions radialWF_;
    ParticlePair entranceParticlePair_;
    MLBWParameters params_;
};

#endif // MLBW_H
