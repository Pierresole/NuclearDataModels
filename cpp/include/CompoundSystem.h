#ifndef COMPOUNDSYSTEM_H
#define COMPOUNDSYSTEM_H

#include <vector>
#include "SpinGroup.h"
#include "ParticlePair.h"
#include <map>
#include <string>
#include <ranges>

class CompoundSystem {
private:
    std::vector<SpinGroup> spinGroups_;
    ParticlePair entranceParticlePair_;

public:
    CompoundSystem(const ParticlePair& entranceParticlePair)
        : entranceParticlePair_(entranceParticlePair)
    {}

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
        std::cout << "Number of spin groups: " << spinGroups_.size() << std::endl;
        for (size_t iSg = 0; iSg < spinGroups_.size(); ++iSg) {
            std::cout << "Spin Group " << iSg << " has " << spinGroups_[iSg].channels().size() << " channels." << std::endl;
        }
    }

    // double computeTotalCrossSection(double E_lab) {
    //     double totalCrossSection = 0.0;
    //     for (const auto& spinGroup : spinGroups_) {
    //         totalCrossSection += spinGroup.computePartialCrossSection(E_lab);
    //     }
    //     return totalCrossSection;
    // }
};

#endif // COMPOUNDSYSTEM_H
