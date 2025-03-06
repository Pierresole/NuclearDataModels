#include "Channel.h"

// void Channel::computeChannelQuantities(double E) {
//     // Compute k^2 using ParticlePair
//     k_c2_ = particlePair_.k2(E);

//     // Compute k
//     k_c_ = std::sqrt(k_c2_);

//     // Compute rho (dimensionless channel radius parameter)
//     double rho = k_c_ * radius_;

//     // Compute P, S, and phi using RadialWaveFunctions
//     radialWaveFunctions_.P_S_Phi(l_, rho);

//     P_c_ = radialWaveFunctions_.getP(l_);
//     S_c_ = radialWaveFunctions_.getS(l_);
//     phi_c_ = std::atan2(radialWaveFunctions_.getSinPhi(l_), radialWaveFunctions_.getCosPhi(l_));
// }

// Compute Shift Factor S_l for l = 0, 1, 2, 3
double Channel::ShiftFactor(double energy, unsigned int l, const ParticlePair& entranceParticlePair) const {

    double rho2 = particlePair_.k2(energy, entranceParticlePair) * std::pow(trueChannelRadius_, 2);
    if (l == 0) return 0.0; // S_0
    if (l == 1) return -1.0 / (1.0 + rho2);
    if (l == 2) return -(18.0 + 3.0 * rho2) / (9.0 + 3.0 * rho2 + std::pow(rho2, 2));
    if (l == 3) return -(675.0 + 90.0 * rho2 + 6.0 * std::pow(rho2, 2)) / 
                        (225.0 + 45.0 * rho2 + 6.0 * std::pow(rho2, 2) + std::pow(rho2, 3));
    throw std::invalid_argument("ShiftFactor only implemented for l = 0, 1, 2, 3.");
}

// Compute Penetration Factor P_l for l = 0, 1, 2, 3
double Channel::PenetrationFactor(double energy, unsigned int l, const ParticlePair& entranceParticlePair) const {
    double rho2 = particlePair_.k2(energy, entranceParticlePair) * std::pow(trueChannelRadius_, 2);
    double rho = std::sqrt(rho2);
    if (l == 0) return rho; // P_0
    if (l == 1) return std::pow(rho, 3) / (1.0 + rho2);
    if (l == 2) return std::pow(rho, 5) / (9.0 + 3.0 * rho2 + std::pow(rho2, 2));
    if (l == 3) return std::pow(rho, 7) / (225.0 + 45.0 * rho2 + 6.0 * std::pow(rho2, 2) + std::pow(rho2, 3));
    throw std::invalid_argument("PenetrationFactor only implemented for l = 0, 1, 2, 3.");
}

// Compute Phase Shift phi_l for l = 0, 1, 2, 3
double Channel::PhaseShift(double energy, unsigned int l, const ParticlePair& entranceParticlePair) const {
    double rho = std::sqrt(particlePair_.k2(energy, entranceParticlePair)) * effectiveChannelRadius_;
    if (l == 0) return rho; // phi_0
    if (l == 1) return rho - std::atan(rho);
    if (l == 2) return rho - std::atan(3.0 * rho / (3.0 - std::pow(rho, 2)));
    if (l == 3) return rho - std::atan(rho * (15.0 - std::pow(rho, 2)) / (15.0 - 6.0 * std::pow(rho, 2)));
    throw std::invalid_argument("PhaseShift only implemented for l = 0, 1, 2, 3.");
}
