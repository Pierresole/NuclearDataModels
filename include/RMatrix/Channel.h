#ifndef CHANNEL_H
#define CHANNEL_H

#include <stdexcept>
#include <tuple>

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

    double k_c_;    // Wave number
    double k_c2_;   // Wave number squared

    // Method to verify channel attributes
    void verifyAttributes() const {
        if (effectiveChannelRadius_ < 0.0) {
            throw std::invalid_argument("Effective channel radius must be positive");
        }
        if (trueChannelRadius_ < 0.0) {
            throw std::invalid_argument("True channel radius must be positive");
        }
        // Add other validations as needed
    }

public:
    Channel(const ParticlePair& particlePair, unsigned int l, double effectiveChannelRadius, double trueChannelRadius, double channelSpin)
        : particlePair_(particlePair), l_(l), effectiveChannelRadius_(effectiveChannelRadius), trueChannelRadius_(trueChannelRadius), spinChannel_(channelSpin)
    {
        verifyAttributes();
    }

    // Return type for channel quantities
    struct ChannelQuantities {
        double P;    // Penetration factor
        double S;    // Shift factor
        double phi;  // Hard-sphere phase shift
        double k2;   // Wave number squared
    };

    // Remove these getters or update them to use computeChannelQuantities
    // double getPenetrationFactor() const { return P_c_; }
    // double getShiftFactor() const { return S_c_; }
    // double getPhaseShift() const { return phi_c_; }
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

    // Updated to return a tuple of values instead of setting attributes
    ChannelQuantities computeChannelQuantities(double E, const ParticlePair& entranceParticlePair) const {
        if(particlePair_.MT() == 102 || particlePair_.MT() == 18) return {1, 0, 0, 0};

        double k2 = particlePair_.k2(E, entranceParticlePair, false);
        double k = std::sqrt(k2);
        
        double rho2 = k2 * std::pow(effectiveChannelRadius_, 2);

        // Use RadialWaveFunctions directly
        RadialWaveFunctions waveFunctions(effectiveChannelRadius_);
        waveFunctions.P_S_Phi(l_, rho2);

        double P = waveFunctions.getP(l_);
        double S = waveFunctions.getS(l_);
        double phi = std::atan2(waveFunctions.getSinPhi(l_), waveFunctions.getCosPhi(l_));
        
        return {P, S, phi, k2};
    }

    double computePenetration(double E, const ParticlePair& entranceParticlePair) const {
        if(particlePair_.MT() == 102 || particlePair_.MT() == 18) return 1;
        // Compute k^2 and k
        double k2 = particlePair_.k2(E, entranceParticlePair, false);
        // double k = std::sqrt(k2);
        // k_c_ = k;    // We still store k for other methods that might need it
        // k_c2_ = k2;  // Same for k2
        
        double rho2 = std::abs(k2) * std::pow(effectiveChannelRadius_, 2);

        // Use RadialWaveFunctions directly
        RadialWaveFunctions waveFunctions(effectiveChannelRadius_);
        return waveFunctions.channelPenetrationAndShift(l_, rho2);
    }

    // Helper method to check if this channel is valid
    bool isValid() const {
        try {
            verifyAttributes();
            return true;
        } catch (const std::invalid_argument&) {
            return false;
        }
    }
};

#endif // CHANNEL_H
