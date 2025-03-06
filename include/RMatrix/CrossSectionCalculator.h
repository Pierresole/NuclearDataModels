#ifndef CROSS_SECTION_CALCULATOR_H
#define CROSS_SECTION_CALCULATOR_H

#include "SpinGroup.h"
#include "RMatrix.h"
#include "LevelMatrix.h"
#include "XMatrix.h"
#include "CollisionMatrix.h"
#include <vector>
#include <string>

class CrossSectionCalculator {
public:
    CrossSectionCalculator(const std::vector<SpinGroup>& spinGroups, double A);

    // Modify the method signature to return vector<double>
    std::vector<double> computeCrossSections(const std::vector<double>& energies, const std::string& reactionChannel);

private:
    std::vector<SpinGroup> spinGroups_;
    double A_;

    // Helper method to check if a channel matches the desired reaction
    bool channelMatches(const Channel& channel, const std::string& reactionChannel) const;
};

#endif // CROSS_SECTION_CALCULATOR_H
