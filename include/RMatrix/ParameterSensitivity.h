#ifndef PARAMETER_SENSITIVITY_H
#define PARAMETER_SENSITIVITY_H

#include "CompoundSystem.h"
#include "SpinGroup.h"
#include <vector>
#include <string>
#include <functional>
#include <map>

// Structure to hold parameter information
struct ParameterInfo {
    enum class Type { RESONANCE_ENERGY, RESONANCE_GAMMA, CHANNEL_RADIUS, PARTICLE_MASS, PARTICLE_SPIN };
    
    Type type;
    size_t spinGroupIndex;
    size_t resonanceIndex;  // For resonance parameters
    size_t channelIndex;    // For channel/gamma parameters
    size_t particlePairIndex; // For particle pair parameters
    std::string description;
    double nominalValue;
    
    ParameterInfo(Type t, size_t sgIdx, size_t resIdx, size_t chIdx, size_t ppIdx, 
                  const std::string& desc, double value)
        : type(t), spinGroupIndex(sgIdx), resonanceIndex(resIdx), 
          channelIndex(chIdx), particlePairIndex(ppIdx), description(desc), nominalValue(value) {}
};

// Structure to define energy groups for multigroup calculations
struct EnergyGroup {
    double lowerBound;
    double upperBound;
    size_t numIntegrationPoints;
    std::string integrationMethod; // "trapezoid", "simpson", "gauss"
    
    EnergyGroup(double lower, double upper, size_t nPoints = 100, const std::string& method = "trapezoid")
        : lowerBound(lower), upperBound(upper), numIntegrationPoints(nPoints), integrationMethod(method) {}
};

class ParameterSensitivity {
private:
    CompoundSystem baseSystem_;
    std::vector<ParameterInfo> parameters_;
    double defaultPerturbation_;
    
public:
    ParameterSensitivity(const CompoundSystem& system, double perturbation = 1e-6)
        : baseSystem_(system), defaultPerturbation_(perturbation) {
        extractParameters();
    }
    
    // Extract all parameters from the compound system
    void extractParameters();
    
    // Compute numerical derivative for a specific parameter
    double computeDerivative(size_t parameterIndex, double energy, 
                           std::function<double(const CompoundSystem&, double)> crossSectionFunc,
                           double perturbation = -1.0);
    
    // Compute derivatives for all parameters
    std::vector<double> computeAllDerivatives(double energy,
                                            std::function<double(const CompoundSystem&, double)> crossSectionFunc,
                                            double perturbation = -1.0);
    
    // Specialized methods for different cross section types
    std::vector<double> totalCrossSectionDerivatives(double energy, double perturbation = -1.0);
    std::vector<double> elasticCrossSectionDerivatives(double energy, double perturbation = -1.0);
    std::vector<double> captureCrossSectionDerivatives(double energy, double perturbation = -1.0);
    std::vector<double> fissionCrossSectionDerivatives(double energy, double perturbation = -1.0);
    
    // New multigroup sensitivity methods
    std::vector<double> multigroupTotalCrossSectionDerivatives(
        const std::vector<EnergyGroup>& energyGroups, 
        double perturbation = -1.0);
    
    std::vector<double> multigroupElasticCrossSectionDerivatives(
        const std::vector<EnergyGroup>& energyGroups, 
        double perturbation = -1.0);
    
    std::vector<double> multigroupCaptureCrossSectionDerivatives(
        const std::vector<EnergyGroup>& energyGroups, 
        double perturbation = -1.0);
    
    std::vector<double> multigroupFissionCrossSectionDerivatives(
        const std::vector<EnergyGroup>& energyGroups, 
        double perturbation = -1.0);
    
    // Generic multigroup derivative computation
    std::vector<double> computeMultigroupDerivatives(
        const std::vector<EnergyGroup>& energyGroups,
        std::function<double(const CompoundSystem&, double)> crossSectionFunc,
        double perturbation = -1.0);
    
    // Multigroup cross section integration (without derivatives)
    std::vector<double> computeMultigroupCrossSections(
        const std::vector<EnergyGroup>& energyGroups,
        std::function<double(const CompoundSystem&, double)> crossSectionFunc) const;
    
    // Enhanced sensitivity matrix computation with multigroup support
    std::map<std::string, std::vector<std::vector<double>>> computeMultigroupSensitivityMatrix(
        const std::vector<EnergyGroup>& energyGroups,
        const std::vector<std::string>& crossSectionTypes = {"total", "elastic", "capture", "fission"});
    
    // Parameter access
    const std::vector<ParameterInfo>& getParameters() const { return parameters_; }
    size_t getParameterCount() const { return parameters_.size(); }
    
    // Create perturbed system
    CompoundSystem createPerturbedSystem(size_t parameterIndex, double perturbation) const;
    
    // Utility methods
    void printParameterInfo() const;
    static std::vector<EnergyGroup> createUniformEnergyGroups(
        double minEnergy, double maxEnergy, size_t numGroups, 
        size_t pointsPerGroup = 100, const std::string& method = "trapezoid");
    
    static std::vector<EnergyGroup> createLogarithmicEnergyGroups(
        double minEnergy, double maxEnergy, size_t numGroups, 
        size_t pointsPerGroup = 100, const std::string& method = "trapezoid");
    
    // Import energy group structure from common formats
    static std::vector<EnergyGroup> fromENDFEnergyGrid(const std::vector<double>& energyBounds);
    static std::vector<EnergyGroup> fromVITAMINJ(const std::string& filename = "");
    
private:
    // Helper methods for parameter modification
    void perturbParameter(CompoundSystem& system, size_t parameterIndex, double perturbation) const;
    SpinGroup createPerturbedSpinGroup(const SpinGroup& original, const ParameterInfo& param, double perturbation) const;
    Channel createPerturbedChannel(const Channel& original, const ParameterInfo& param, double perturbation) const;
    ParticlePair createPerturbedParticlePair(const ParticlePair& original, const ParameterInfo& param, double perturbation) const;
    Resonance createPerturbedResonance(const Resonance& original, const ParameterInfo& param, double perturbation) const;
    
    // Integration methods
    double integrateOverEnergyGroup(
        const EnergyGroup& group,
        std::function<double(double)> func) const;
    
    double trapezoidalIntegration(
        double lowerBound, double upperBound, size_t numPoints,
        std::function<double(double)> func) const;
    
    double simpsonIntegration(
        double lowerBound, double upperBound, size_t numPoints,
        std::function<double(double)> func) const;
    
    double gaussianQuadrature(
        double lowerBound, double upperBound, size_t numPoints,
        std::function<double(double)> func) const;
    
    // Helper for creating energy point arrays
    std::vector<double> createEnergyPoints(const EnergyGroup& group) const;
};

#endif // PARAMETER_SENSITIVITY_H
