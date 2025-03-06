#ifndef CHANNEL_H
#define CHANNEL_H

#include <stdexcept>

#include "ParticlePair.h"
#include "RadialWaveFunctions.h"

class Channel {
    
private:
    ParticlePair particlePair_;
    unsigned int l_;
    double trueChannelRadius_;
    double effectiveChannelRadius_;
    double spinChannel_; 
    //RadialWaveFunctions radialWaveFunctions_;

    double P_c_;    // Penetration factor
    double S_c_;    // Shift factor
    double phi_c_;  // Hard-sphere phase shift
    double k_c_;    // Wave number
    double k_c2_;   // Wave number squared

public:
    Channel(const ParticlePair& particlePair, unsigned int l, double effectiveChannelRadius, double trueChannelRadius, double channelSpin)
        : particlePair_(particlePair), l_(l), effectiveChannelRadius_(effectiveChannelRadius), trueChannelRadius_(trueChannelRadius), spinChannel_(channelSpin)
    {}

    // void computeChannelQuantities(double E);

    double getPenetrationFactor() const { return P_c_; }
    double getShiftFactor() const { return S_c_; }
    double getPhaseShift() const { return phi_c_; }
    double getWaveNumber() const { return k_c_; }
    double getWaveNumberSquared() const { return k_c2_; }
    double mass1() const { return particlePair_.mass1(); }
    double mass2() const { return particlePair_.mass2(); }
    unsigned int L() const { return l_; }

    const ParticlePair& getParticlePair() const { return particlePair_; }

    double waveNumber() const { return k_c_; }
    double PhaseShift(double energy, unsigned int l, const ParticlePair& entranceParticlePair) const ;
    double PenetrationFactor(double energy, unsigned int l, const ParticlePair& entranceParticlePair) const ;
    double ShiftFactor(double energy, unsigned int l, const ParticlePair& entranceParticlePair) const ;

    void computeChannelQuantities(double E, const ParticlePair& entranceParticlePair) {
        // Compute k^2 and k
        k_c2_ = particlePair_.k2(E, entranceParticlePair, false);
        k_c_ = std::sqrt(k_c2_);
        double rho = k_c_ * effectiveChannelRadius_;

        // Use RadialWaveFunctions directly
        RadialWaveFunctions waveFunctions(effectiveChannelRadius_);
        waveFunctions.P_S_Phi(l_, rho);

        P_c_ = waveFunctions.getP(l_);
        S_c_ = waveFunctions.getS(l_);
        phi_c_ = std::atan2(waveFunctions.getSinPhi(l_), waveFunctions.getCosPhi(l_));
    }
};

#endif // CHANNEL_H
