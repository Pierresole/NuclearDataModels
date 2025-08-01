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
            sigma_total += sg.crossSection(E, entranceParticlePair_);
        }
        return sigma_total;
    }
    
    // Elastic cross section summed over all spin groups
    double elasticCrossSection(double E) const {
        double sigma_elastic = 0.0;
        for (const auto& sg : spinGroups_) {
            sigma_elastic += sg.elasticCrossSection(E, entranceParticlePair_);
        }
        return sigma_elastic;
    }
    
    // Capture cross section summed over all spin groups
    double captureCrossSection(double E) const {
        double sigma_capture = 0.0;
        for (const auto& sg : spinGroups_) {
            sigma_capture += sg.captureCrossSection(E, entranceParticlePair_);
        }
        return sigma_capture;
    }
    
    // Fission cross section summed over all spin groups
    double fissionCrossSection(double E) const {
        double sigma_fission = 0.0;
        for (const auto& sg : spinGroups_) {
            sigma_fission += sg.fissionCrossSection(E, entranceParticlePair_);
        }
        return sigma_fission;
    }
    
    // Total cross section summed over all spin groups
    double totalCrossSection(double E) const {
        double sigma_total = 0.0;
        for (const auto& sg : spinGroups_) {
            sigma_total += sg.totalCrossSection(E, entranceParticlePair_);
        }
        return sigma_total;
    }
    
    // Get partial cross section for a specific spin group
    double spinGroupElasticCrossSection(size_t spinGroupIndex, double E) const {
        if (spinGroupIndex >= spinGroups_.size()) {
            throw std::out_of_range("Spin group index out of range");
        }
        const auto& sg = spinGroups_[spinGroupIndex];
        return sg.elasticCrossSection(E, entranceParticlePair_);
    }
    
    double spinGroupCaptureCrossSection(size_t spinGroupIndex, double E) const {
        if (spinGroupIndex >= spinGroups_.size()) {
            throw std::out_of_range("Spin group index out of range");
        }
        const auto& sg = spinGroups_[spinGroupIndex];
        return sg.captureCrossSection(E, entranceParticlePair_);
    }
    
    double spinGroupFissionCrossSection(size_t spinGroupIndex, double E) const {
        if (spinGroupIndex >= spinGroups_.size()) {
            throw std::out_of_range("Spin group index out of range");
        }
        const auto& sg = spinGroups_[spinGroupIndex];
        return sg.fissionCrossSection(E, entranceParticlePair_);
    }
    
    double spinGroupTotalCrossSection(size_t spinGroupIndex, double E) const {
        if (spinGroupIndex >= spinGroups_.size()) {
            throw std::out_of_range("Spin group index out of range");
        }
        const auto& sg = spinGroups_[spinGroupIndex];
        return sg.totalCrossSection(E, entranceParticlePair_);
    }
};

#endif // COMPOUNDSYSTEM_H
