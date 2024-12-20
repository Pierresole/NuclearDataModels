#ifndef PARTICLE_PAIR_H
#define PARTICLE_PAIR_H

#include <string>

class ParticlePair {
private:
    double massParticleA_; // Mass of neutron or incident particle
    double massParticleB_;
    double spinParticleA_;
    double spinParticleB_;
    double QI_;
    int parityParticleA_;
    int parityParticleB_;
    bool entrance_;
    int MT_;
    std::string reactionID_;

public:
    ParticlePair(double massA, double massB, double spinA, double spinB, double QI,
                 int parityA, int parityB, int MT)
        : massParticleA_(massA), massParticleB_(massB),
          spinParticleA_(spinA), spinParticleB_(spinB),
          QI_(QI), parityParticleA_(parityA), parityParticleB_(parityB), MT_(MT) {}

    ParticlePair(double massB, const std::string& reactionID) : massParticleB_(massB), reactionID_(reactionID) {};
    ParticlePair(double massB) : massParticleB_(massB) {};

    // Method to calculate wave number squared (k^2)
    double k2(double E, const ParticlePair& entranceParticlePair, bool isKinematicRelativistic = false) const {
        double constantMomentum2WaveNumber2 = 2.392253439955726e-6;
        return 2 * reducedMass() * (entranceParticlePair.mass2() / (entranceParticlePair.mass2() + entranceParticlePair.mass1() ) * E + QI_) * constantMomentum2WaveNumber2;
    }

    // Getter for reaction ID
    const std::string& getReactionID() const { return reactionID_; }

    double reducedMass() const { return massParticleA_ * massParticleB_ / (massParticleA_ + massParticleB_); }

    double mass1() const { return massParticleA_; }
    double mass2() const { return massParticleB_; }
    double spin1() const { return spinParticleA_; }
    double spin2() const { return spinParticleB_; }
    int MT() const { return MT_; }
};

#endif // PARTICLE_PAIR_H
