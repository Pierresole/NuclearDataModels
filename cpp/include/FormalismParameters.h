// FormalismParameters.h.h
#ifndef FORMALISMPARAMETERS_H
#define FORMALISMPARAMETERS_H

#include <vector>
#include <variant>
#include <optional>

struct RPlist {
    std::vector<double> ER;
    std::vector<double> AJ;
    std::vector<double> GN;
    std::vector<double> GG;
    std::vector<double> GFA;
    std::vector<double> GFB;
};

struct Llist {
    double APL;
    int L;
    RPlist resonanceParameters;  // Nested structure containing resonance parameters
};

struct JPiList {
    double J;
    double PI;
    // Add additional R-Matrix parameters as needed
};

struct RMatrixParameters {
    std::vector<JPiList> JPiGroups;  // R-Matrix specific parameter group
};

struct MLBWParameters {
    double SPI;
    double AP;
    double AWRI;
    std::vector<Llist> lValues;  // List of L values for MLBW
};

// Variant to hold either MLBW or R-Matrix parameters
using ResonanceParameters = std::variant<MLBWParameters, RMatrixParameters>;

struct FormalismParameters {
    double SPI;
    double AP;
    int LAD;
    ResonanceParameters parameters;  // Can hold MLBW or R-Matrix parameters
};

#endif // FORMALISMPARAMETERS_H
