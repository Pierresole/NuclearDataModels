#include "ParameterSensitivity.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>

void ParameterSensitivity::extractParameters() {
    parameters_.clear();
    
    auto spinGroups = baseSystem_.spinGroups();
    size_t sgIndex = 0;
    
    for (const auto& sg : spinGroups) {
        // Extract resonance parameters
        const auto& resonances = sg.getResonances();
        for (size_t resIdx = 0; resIdx < resonances.size(); ++resIdx) {
            const auto& res = resonances[resIdx];
            
            // Resonance energy
            parameters_.emplace_back(
                ParameterInfo::Type::RESONANCE_ENERGY,
                sgIndex, resIdx, 0, 0,
                "SG" + std::to_string(sgIndex) + "_RES" + std::to_string(resIdx) + "_Energy",
                res.getEnergy()
            );
            
            // Resonance gammas (one for each channel)
            const auto& gammas = res.getGamma();
            for (size_t chIdx = 0; chIdx < gammas.size(); ++chIdx) {
                parameters_.emplace_back(
                    ParameterInfo::Type::RESONANCE_GAMMA,
                    sgIndex, resIdx, chIdx, 0,
                    "SG" + std::to_string(sgIndex) + "_RES" + std::to_string(resIdx) + "_GAMMA" + std::to_string(chIdx),
                    gammas[chIdx]
                );
            }
        }
        
        // Extract channel parameters (radii)
        const auto& channels = sg.channels();
        for (size_t chIdx = 0; chIdx < channels.size(); ++chIdx) {
            // Note: Channel radius extraction would need additional methods in Channel class
            // For now, we'll focus on resonance parameters
        }
        
        sgIndex++;
    }
    
    std::cout << "Extracted " << parameters_.size() << " parameters for sensitivity analysis." << std::endl;
}

double ParameterSensitivity::computeDerivative(size_t parameterIndex, double energy,
                                             std::function<double(const CompoundSystem&, double)> crossSectionFunc,
                                             double perturbation) {
    if (perturbation < 0) perturbation = defaultPerturbation_;
    
    if (parameterIndex >= parameters_.size()) {
        throw std::out_of_range("Parameter index out of range");
    }
    
    const auto& param = parameters_[parameterIndex];
    
    // Use relative perturbation for better numerical stability
    double actualPerturbation = perturbation * std::abs(param.nominalValue);
    if (actualPerturbation == 0.0) actualPerturbation = perturbation; // Fallback for zero parameters
    
    // Forward difference approximation
    CompoundSystem perturbedSystem = createPerturbedSystem(parameterIndex, actualPerturbation);
    
    double originalValue = crossSectionFunc(baseSystem_, energy);
    double perturbedValue = crossSectionFunc(perturbedSystem, energy);
    
    return (perturbedValue - originalValue) / actualPerturbation;
}

std::vector<double> ParameterSensitivity::computeAllDerivatives(double energy,
                                                              std::function<double(const CompoundSystem&, double)> crossSectionFunc,
                                                              double perturbation) {
    std::vector<double> derivatives;
    derivatives.reserve(parameters_.size());
    
    for (size_t i = 0; i < parameters_.size(); ++i) {
        derivatives.push_back(computeDerivative(i, energy, crossSectionFunc, perturbation));
    }
    
    return derivatives;
}

std::vector<double> ParameterSensitivity::totalCrossSectionDerivatives(double energy, double perturbation) {
    return computeAllDerivatives(energy, 
        [](const CompoundSystem& sys, double E) { return sys.totalCrossSection(E); }, 
        perturbation);
}

std::vector<double> ParameterSensitivity::elasticCrossSectionDerivatives(double energy, double perturbation) {
    return computeAllDerivatives(energy,
        [](const CompoundSystem& sys, double E) { return sys.elasticCrossSection(E); },
        perturbation);
}

std::vector<double> ParameterSensitivity::captureCrossSectionDerivatives(double energy, double perturbation) {
    return computeAllDerivatives(energy,
        [](const CompoundSystem& sys, double E) { return sys.captureCrossSection(E); },
        perturbation);
}

std::vector<double> ParameterSensitivity::fissionCrossSectionDerivatives(double energy, double perturbation) {
    return computeAllDerivatives(energy,
        [](const CompoundSystem& sys, double E) { return sys.fissionCrossSection(E); },
        perturbation);
}

CompoundSystem ParameterSensitivity::createPerturbedSystem(size_t parameterIndex, double perturbation) const {
    if (parameterIndex >= parameters_.size()) {
        throw std::out_of_range("Parameter index out of range");
    }
    
    // Create a copy of the base system
    CompoundSystem perturbedSystem = baseSystem_;
    
    // Apply perturbation
    perturbParameter(perturbedSystem, parameterIndex, perturbation);
    
    return perturbedSystem;
}

void ParameterSensitivity::perturbParameter(CompoundSystem& system, size_t parameterIndex, double perturbation) const {
    const auto& param = parameters_[parameterIndex];
    
    // This is a simplified approach - in practice, you'd need to reconstruct the system
    // with perturbed parameters. For now, we'll create a new system from scratch.
    
    // Get the original spin groups
    auto originalSpinGroups = baseSystem_.spinGroups();
    std::vector<SpinGroup> newSpinGroups;
    
    size_t sgIdx = 0;
    for (const auto& sg : originalSpinGroups) {
        if (sgIdx == param.spinGroupIndex) {
            // This spin group contains the parameter to perturb
            newSpinGroups.push_back(createPerturbedSpinGroup(sg, param, perturbation));
        } else {
            // Copy unchanged
            newSpinGroups.push_back(sg);
        }
        sgIdx++;
    }
    
    // Reconstruct the system
    system = CompoundSystem(baseSystem_.entranceParticlePair(), newSpinGroups);
}

SpinGroup ParameterSensitivity::createPerturbedSpinGroup(const SpinGroup& original, const ParameterInfo& param, double perturbation) const {
    if (param.type == ParameterInfo::Type::RESONANCE_ENERGY || 
        param.type == ParameterInfo::Type::RESONANCE_GAMMA) {
        
        // Create new resonances with perturbation
        std::vector<Resonance> newResonances;
        const auto& originalResonances = original.getResonances();
        
        for (size_t resIdx = 0; resIdx < originalResonances.size(); ++resIdx) {
            if (resIdx == param.resonanceIndex) {
                newResonances.push_back(createPerturbedResonance(originalResonances[resIdx], param, perturbation));
            } else {
                newResonances.push_back(originalResonances[resIdx]);
            }
        }
        
        // Create new spin group with perturbed resonances
        return SpinGroup(original.getJ(), original.getPJ(), original.channels(), newResonances);
    }
    
    // For other parameter types, return original (to be implemented)
    return original;
}

Resonance ParameterSensitivity::createPerturbedResonance(const Resonance& original, const ParameterInfo& param, double perturbation) const {
    if (param.type == ParameterInfo::Type::RESONANCE_ENERGY) {
        return Resonance(original.getEnergy() + perturbation, original.getGamma());
    } else if (param.type == ParameterInfo::Type::RESONANCE_GAMMA) {
        auto newGammas = original.getGamma();
        newGammas[param.channelIndex] += perturbation;
        return Resonance(original.getEnergy(), newGammas);
    }
    
    return original;
}

void ParameterSensitivity::printParameterInfo() const {
    std::cout << std::setw(5) << "Index" 
              << std::setw(20) << "Description" 
              << std::setw(15) << "Type"
              << std::setw(15) << "Nominal Value" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
    
    for (size_t i = 0; i < parameters_.size(); ++i) {
        const auto& param = parameters_[i];
        std::string typeStr;
        switch (param.type) {
            case ParameterInfo::Type::RESONANCE_ENERGY: typeStr = "RES_ENERGY"; break;
            case ParameterInfo::Type::RESONANCE_GAMMA: typeStr = "RES_GAMMA"; break;
            case ParameterInfo::Type::CHANNEL_RADIUS: typeStr = "CH_RADIUS"; break;
            case ParameterInfo::Type::PARTICLE_MASS: typeStr = "PART_MASS"; break;
            case ParameterInfo::Type::PARTICLE_SPIN: typeStr = "PART_SPIN"; break;
        }
        
        std::cout << std::setw(5) << i
                  << std::setw(20) << param.description
                  << std::setw(15) << typeStr
                  << std::setw(15) << param.nominalValue << std::endl;
    }
}

std::map<std::string, std::vector<double>> ParameterSensitivity::computeSensitivityMatrix(
    const std::vector<double>& energies, 
    const std::vector<std::string>& crossSectionTypes) {
    
    std::map<std::string, std::vector<double>> results;
    
    for (const auto& type : crossSectionTypes) {
        std::vector<double> sensitivities;
        sensitivities.reserve(energies.size() * parameters_.size());
        
        for (double energy : energies) {
            std::vector<double> derivatives;
            
            if (type == "total") {
                derivatives = totalCrossSectionDerivatives(energy);
            } else if (type == "elastic") {
                derivatives = elasticCrossSectionDerivatives(energy);
            } else if (type == "capture") {
                derivatives = captureCrossSectionDerivatives(energy);
            } else if (type == "fission") {
                derivatives = fissionCrossSectionDerivatives(energy);
            }
            
            sensitivities.insert(sensitivities.end(), derivatives.begin(), derivatives.end());
        }
        
        results[type] = sensitivities;
    }
    
    return results;
}

std::vector<double> ParameterSensitivity::multigroupTotalCrossSectionDerivatives(
    const std::vector<EnergyGroup>& energyGroups, double perturbation) {
    return computeMultigroupDerivatives(energyGroups,
        [](const CompoundSystem& sys, double E) { return sys.totalCrossSection(E); },
        perturbation);
}

std::vector<double> ParameterSensitivity::multigroupElasticCrossSectionDerivatives(
    const std::vector<EnergyGroup>& energyGroups, double perturbation) {
    return computeMultigroupDerivatives(energyGroups,
        [](const CompoundSystem& sys, double E) { return sys.elasticCrossSection(E); },
        perturbation);
}

std::vector<double> ParameterSensitivity::multigroupCaptureCrossSectionDerivatives(
    const std::vector<EnergyGroup>& energyGroups, double perturbation) {
    return computeMultigroupDerivatives(energyGroups,
        [](const CompoundSystem& sys, double E) { return sys.captureCrossSection(E); },
        perturbation);
}

std::vector<double> ParameterSensitivity::multigroupFissionCrossSectionDerivatives(
    const std::vector<EnergyGroup>& energyGroups, double perturbation) {
    return computeMultigroupDerivatives(energyGroups,
        [](const CompoundSystem& sys, double E) { return sys.fissionCrossSection(E); },
        perturbation);
}

std::vector<double> ParameterSensitivity::computeMultigroupDerivatives(
    const std::vector<EnergyGroup>& energyGroups,
    std::function<double(const CompoundSystem&, double)> crossSectionFunc,
    double perturbation) {
    
    if (perturbation < 0) perturbation = defaultPerturbation_;
    
    size_t numGroups = energyGroups.size();
    size_t numParams = parameters_.size();
    std::vector<double> derivatives(numGroups * numParams, 0.0);
    
    // Compute multigroup cross sections for base system
    std::vector<double> baseMGXS = computeMultigroupCrossSections(energyGroups, crossSectionFunc);
    
    // For each parameter
    for (size_t paramIdx = 0; paramIdx < numParams; ++paramIdx) {
        const auto& param = parameters_[paramIdx];
        
        // Use relative perturbation for better numerical stability
        double actualPerturbation = perturbation * std::abs(param.nominalValue);
        if (actualPerturbation == 0.0) actualPerturbation = perturbation;
        
        // Create perturbed system
        CompoundSystem perturbedSystem = createPerturbedSystem(paramIdx, actualPerturbation);
        
        // Compute multigroup cross sections for perturbed system
        std::vector<double> perturbedMGXS = computeMultigroupCrossSections(energyGroups, 
            [&perturbedSystem, crossSectionFunc](double E) { return crossSectionFunc(perturbedSystem, E); });
        
        // Compute finite difference derivatives for each group
        for (size_t groupIdx = 0; groupIdx < numGroups; ++groupIdx) {
            size_t linearIdx = groupIdx * numParams + paramIdx;
            derivatives[linearIdx] = (perturbedMGXS[groupIdx] - baseMGXS[groupIdx]) / actualPerturbation;
        }
        
        // Progress indicator for large calculations
        if (numParams > 50 && (paramIdx % 10 == 0)) {
            std::cout << "Computed derivatives for parameter " << paramIdx << "/" << numParams << std::endl;
        }
    }
    
    return derivatives;
}

std::vector<double> ParameterSensitivity::computeMultigroupCrossSections(
    const std::vector<EnergyGroup>& energyGroups,
    std::function<double(const CompoundSystem&, double)> crossSectionFunc) const {
    
    std::vector<double> mgCrossSections;
    mgCrossSections.reserve(energyGroups.size());
    
    for (const auto& group : energyGroups) {
        // Define the function to integrate over this energy group
        auto integrand = [this, crossSectionFunc](double E) {
            return crossSectionFunc(baseSystem_, E);
        };
        
        // Integrate over the energy group
        double integral = integrateOverEnergyGroup(group, integrand);
        
        // Normalize by group width to get average cross section
        double groupWidth = group.upperBound - group.lowerBound;
        mgCrossSections.push_back(integral / groupWidth);
    }
    
    return mgCrossSections;
}

std::map<std::string, std::vector<std::vector<double>>> ParameterSensitivity::computeMultigroupSensitivityMatrix(
    const std::vector<EnergyGroup>& energyGroups,
    const std::vector<std::string>& crossSectionTypes) {
    
    std::map<std::string, std::vector<std::vector<double>>> results;
    
    size_t numGroups = energyGroups.size();
    size_t numParams = parameters_.size();
    
    for (const auto& type : crossSectionTypes) {
        std::vector<double> flatDerivatives;
        
        if (type == "total") {
            flatDerivatives = multigroupTotalCrossSectionDerivatives(energyGroups);
        } else if (type == "elastic") {
            flatDerivatives = multigroupElasticCrossSectionDerivatives(energyGroups);
        } else if (type == "capture") {
            flatDerivatives = multigroupCaptureCrossSectionDerivatives(energyGroups);
        } else if (type == "fission") {
            flatDerivatives = multigroupFissionCrossSectionDerivatives(energyGroups);
        }
        
        // Reshape flat array into 2D matrix: [group][parameter]
        std::vector<std::vector<double>> matrix(numGroups, std::vector<double>(numParams));
        for (size_t g = 0; g < numGroups; ++g) {
            for (size_t p = 0; p < numParams; ++p) {
                matrix[g][p] = flatDerivatives[g * numParams + p];
            }
        }
        
        results[type] = matrix;
    }
    
    return results;
}

double ParameterSensitivity::integrateOverEnergyGroup(
    const EnergyGroup& group, std::function<double(double)> func) const {
    
    if (group.integrationMethod == "trapezoid") {
        return trapezoidalIntegration(group.lowerBound, group.upperBound, 
                                    group.numIntegrationPoints, func);
    } else if (group.integrationMethod == "simpson") {
        return simpsonIntegration(group.lowerBound, group.upperBound, 
                                group.numIntegrationPoints, func);
    } else if (group.integrationMethod == "gauss") {
        return gaussianQuadrature(group.lowerBound, group.upperBound, 
                                group.numIntegrationPoints, func);
    } else {
        // Default to trapezoidal
        return trapezoidalIntegration(group.lowerBound, group.upperBound, 
                                    group.numIntegrationPoints, func);
    }
}

double ParameterSensitivity::trapezoidalIntegration(
    double lowerBound, double upperBound, size_t numPoints,
    std::function<double(double)> func) const {
    
    if (numPoints < 2) numPoints = 2;
    
    double h = (upperBound - lowerBound) / (numPoints - 1);
    double sum = 0.5 * (func(lowerBound) + func(upperBound));
    
    for (size_t i = 1; i < numPoints - 1; ++i) {
        double x = lowerBound + i * h;
        sum += func(x);
    }
    
    return sum * h;
}

double ParameterSensitivity::simpsonIntegration(
    double lowerBound, double upperBound, size_t numPoints,
    std::function<double(double)> func) const {
    
    // Ensure even number of intervals for Simpson's rule
    if (numPoints % 2 == 0) numPoints++;
    if (numPoints < 3) numPoints = 3;
    
    double h = (upperBound - lowerBound) / (numPoints - 1);
    double sum = func(lowerBound) + func(upperBound);
    
    for (size_t i = 1; i < numPoints - 1; i += 2) {
        double x = lowerBound + i * h;
        sum += 4.0 * func(x);
    }
    
    for (size_t i = 2; i < numPoints - 1; i += 2) {
        double x = lowerBound + i * h;
        sum += 2.0 * func(x);
    }
    
    return sum * h / 3.0;
}

double ParameterSensitivity::gaussianQuadrature(
    double lowerBound, double upperBound, size_t numPoints,
    std::function<double(double)> func) const {
    
    // Simple Gaussian quadrature implementation for common point counts
    // For a more complete implementation, you'd want to use a library like Boost
    
    if (numPoints == 1) {
        return (upperBound - lowerBound) * func((lowerBound + upperBound) / 2.0);
    } else if (numPoints == 2) {
        double a = (upperBound - lowerBound) / 2.0;
        double b = (upperBound + lowerBound) / 2.0;
        double x1 = b - a / std::sqrt(3.0);
        double x2 = b + a / std::sqrt(3.0);
        return a * (func(x1) + func(x2));
    } else {
        // For higher order, fall back to trapezoidal
        return trapezoidalIntegration(lowerBound, upperBound, numPoints, func);
    }
}

std::vector<EnergyGroup> ParameterSensitivity::createUniformEnergyGroups(
    double minEnergy, double maxEnergy, size_t numGroups, 
    size_t pointsPerGroup, const std::string& method) {
    
    std::vector<EnergyGroup> groups;
    groups.reserve(numGroups);
    
    double groupWidth = (maxEnergy - minEnergy) / numGroups;
    
    for (size_t i = 0; i < numGroups; ++i) {
        double lower = minEnergy + i * groupWidth;
        double upper = minEnergy + (i + 1) * groupWidth;
        groups.emplace_back(lower, upper, pointsPerGroup, method);
    }
    
    return groups;
}

std::vector<EnergyGroup> ParameterSensitivity::createLogarithmicEnergyGroups(
    double minEnergy, double maxEnergy, size_t numGroups, 
    size_t pointsPerGroup, const std::string& method) {
    
    std::vector<EnergyGroup> groups;
    groups.reserve(numGroups);
    
    double logMin = std::log(minEnergy);
    double logMax = std::log(maxEnergy);
    double logWidth = (logMax - logMin) / numGroups;
    
    for (size_t i = 0; i < numGroups; ++i) {
        double logLower = logMin + i * logWidth;
        double logUpper = logMin + (i + 1) * logWidth;
        double lower = std::exp(logLower);
        double upper = std::exp(logUpper);
        groups.emplace_back(lower, upper, pointsPerGroup, method);
    }
    
    return groups;
}

std::vector<EnergyGroup> ParameterSensitivity::fromENDFEnergyGrid(const std::vector<double>& energyBounds) {
    std::vector<EnergyGroup> groups;
    groups.reserve(energyBounds.size() - 1);
    
    for (size_t i = 0; i < energyBounds.size() - 1; ++i) {
        groups.emplace_back(energyBounds[i], energyBounds[i + 1]);
    }
    
    return groups;
}

std::vector<double> ParameterSensitivity::createEnergyPoints(const EnergyGroup& group) const {
    std::vector<double> points;
    points.reserve(group.numIntegrationPoints);
    
    double step = (group.upperBound - group.lowerBound) / (group.numIntegrationPoints - 1);
    
    for (size_t i = 0; i < group.numIntegrationPoints; ++i) {
        points.push_back(group.lowerBound + i * step);
    }
    
    return points;
}

// Additional helper methods would be implemented here...
