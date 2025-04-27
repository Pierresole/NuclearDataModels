#ifndef COMPOUNDSYSTEM_H
#define COMPOUNDSYSTEM_H

#include <vector>
#include "SpinGroup.h"
#include "ParticlePair.h"
#include <map>
#include <iostream>
#include <string>
#include <ranges>

class CompoundSystem {
private:
    std::vector<SpinGroup> spinGroups_;
    ParticlePair entranceParticlePair_;

public:
    CompoundSystem(const ParticlePair& entranceParticlePair, const std::vector<SpinGroup>& spinGroups)
        : entranceParticlePair_(entranceParticlePair), spinGroups_(spinGroups)
    {}


    void addSpinGroup(SpinGroup& spinGroup) {
        spinGroups_.push_back(spinGroup);
    }


    auto spinGroups() const { 
        return std::ranges::views::all( this->spinGroups_ ); 
    }
    
    const ParticlePair& entranceParticlePair() const { return entranceParticlePair_; }
    
    const SpinGroup& getSpinGroup(size_t index) const {
        return spinGroups_.at(index);
    }

    void printSpinGroupInfo() const {
        for (size_t iSg = 0; iSg < spinGroups_.size(); ++iSg) {
            std::cout << "Spin Group " << iSg << "/" << spinGroups_.size() << " (" << spinGroups_[iSg].getJ() << ", "<< spinGroups_[iSg].getPJ() << "):"
                      << " has " << spinGroups_[iSg].channels().size() << " channels." << std::endl;
            for(size_t iRes = 0; iRes < spinGroups_[iSg].getResonances().size(); iRes++){
                std::cout << "Resonance " << spinGroups_[iSg].getResonances()[iRes].getEnergy() << " ";
                for(size_t iCH = 0; iCH < spinGroups_[iSg].channels().size(); iCH++){
                    std::cout << spinGroups_[iSg].getResonances()[iRes].getGamma()[iCH] << " ";
                }
                std::cout << std::endl;
            }
        }
    }
    
    double crossSection(double E) const {
        double sigma_total = 0.0;
        for (const auto& sg : spinGroups_) {
            sigma_total += (2 * sg.getJ() + 1) / 
                        ((2 * entranceParticlePair_.spin1() + 1)*(2 * entranceParticlePair_.spin2() + 1)) * 
                        sg.crossSection(E, entranceParticlePair_);
        }
        return sigma_total;
    }
};

#endif // COMPOUNDSYSTEM_H
