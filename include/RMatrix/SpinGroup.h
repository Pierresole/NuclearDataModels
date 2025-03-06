#ifndef SPIN_GROUP_H
#define SPIN_GROUP_H

#include "Channel.h"
#include "ParticlePair.h"
#include <vector>
#include <memory>

class Resonance {
private:
    double energy_;
    std::vector<double> gamma_; // Reduced-width amplitudes for each channel
public:
    Resonance(double energy, const std::vector<double>& gamma)
        : energy_(energy), gamma_(gamma) {}

    double getEnergy() const { return energy_; }
    const std::vector<double>& getGamma() const { return gamma_; }

};

/* SpinGroup class
*  A group of channels with same total angular momentum and parity is called a spin group.
*  The resonances are defined by their coupling to all channels in the spin group via the reduced widths.
*
*/
class SpinGroup {
    
private:
    double J_;
    int PJ_;
    std::vector<Channel> channels_;
    std::vector<Resonance> resonances_;

public:
    SpinGroup(double J, int PJ)
        : J_(J), PJ_(PJ) {}

    SpinGroup(double J, int PJ, const std::vector<Channel>& channels, const std::vector<Resonance>& resonances)
        : J_(J), PJ_(PJ), channels_(channels), resonances_(resonances) {}

    void addChannel(const Channel& channel) { channels_.push_back(channel); }
    void addResonance(const Resonance& resonance) { resonances_.push_back(resonance); }

    const std::vector<Channel>& channels() const { return channels_; }
    const std::vector<Resonance>& getResonances() const { return resonances_; }
    double getJ() const { return J_; }
    int getPJ() const { return PJ_; }

    // double computePartialCrossSection(double E_com) const {
    //     double partialCrossSection = 0.0;

    //     // For each channel in the spin group
    //     for (const auto& channel : channels_) {
    //         // Compute channel quantities at E_com
    //         //channel.computeChannelQuantities(E_com);

    //         // Compute partial R-matrix, U-matrix, etc.
    //         // ...computations...

    //         // Sum up the partial cross sections
    //         // partialCrossSection += ...;
    //     }

    //     return partialCrossSection;
    // }
};

#endif // SPIN_GROUP_H
