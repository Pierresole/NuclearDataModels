// RadialWaveFunctions.h
#ifndef RADIAL_WAVE_FUNCTIONS_H
#define RADIAL_WAVE_FUNCTIONS_H

#include <vector>
#include <cmath>

class RadialWaveFunctions {
public:
    RadialWaveFunctions(double AP) : AP(AP), rho2_(0.0) {}

    void P_S_Phi(unsigned int iL, double dRho2, double dEta = 0.0) const;

    double getP(unsigned int l) const { return P_.at(l); }
    double getS(unsigned int l) const { return S_.at(l); }
    double getCosPhi(unsigned int l) const { return cos_phi_.at(l); }
    double getSinPhi(unsigned int l) const { return sin_phi_.at(l); }
    double getCos2Phi(unsigned int l) const { return cos_2phi_.at(l); }
    double getSin2Phi(unsigned int l) const { return sin_2phi_.at(l); }

private:
    double AP;
    mutable double rho2_;

    mutable std::vector<double> P_;
    mutable std::vector<double> S_;
    mutable std::vector<double> cos_phi_;
    mutable std::vector<double> sin_phi_;
    mutable std::vector<double> cos_2phi_;
    mutable std::vector<double> sin_2phi_;

    void computePenetration(unsigned int l, double dRho2) const;
    void computeShift(unsigned int l, double dRho2) const;
    void computePhase(unsigned int l) const;
};

#endif // RADIAL_WAVE_FUNCTIONS_H
