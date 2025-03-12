// RadialWaveFunctions.cpp
#include "RadialWaveFunctions.h"

void RadialWaveFunctions::P_S_Phi(unsigned int iL, double dRho2, double dEta) const {
    if (dRho2 != rho2_) {
        P_.clear();
        S_.clear();
        cos_phi_.clear();
        sin_phi_.clear();
        cos_2phi_.clear();
        sin_2phi_.clear();
        rho2_ = dRho2;
    }

    if (dRho2 >= 0) {
        while (P_.size() < iL + 1) {
            unsigned int l = P_.size();

            switch (l) {
                case 0:
                    P_.push_back(std::sqrt(dRho2));
                    S_.push_back(0.0);
                    cos_phi_.push_back(std::cos(std::sqrt(dRho2)));
                    sin_phi_.push_back(std::sin(std::sqrt(dRho2)));
                    cos_2phi_.push_back(std::cos(2.0 * std::sqrt(dRho2)));
                    sin_2phi_.push_back(std::sin(2.0 * std::sqrt(dRho2)));
                    break;

                default:
                    double dP = P_.back();
                    double dS = l - S_.back();
                    double dY = dP / dS;
                    double dX = 1.0 - std::pow(dY, 2);
                    double dN = 1.0 / (1.0 + std::pow(dY, 2));
                    double dZ = std::sqrt(dN);
                    double dC = dRho2 / (std::pow(dS, 2) + std::pow(dP, 2));
                    double dCos = cos_phi_.back();
                    double dSin = sin_phi_.back();
                    double dCos_2 = cos_2phi_.back();
                    double dSin_2 = sin_2phi_.back();

                    P_.push_back(dC * dP);
                    S_.push_back(dC * dS - l);

                    cos_phi_.push_back((dCos + dY * dSin) * dZ);
                    sin_phi_.push_back((dSin - dY * dCos) * dZ);
                    cos_2phi_.push_back((dX * dCos_2 + 2 * dY * dSin_2) * dN);
                    sin_2phi_.push_back((dX * dSin_2 - 2 * dY * dCos_2) * dN);

                    break;
            }
        }
    } else {
        while (P_.size() < iL + 1) {
            unsigned int l = P_.size();

            switch (l) {
                case 0:
                    P_.push_back(0.0);
                    S_.push_back(-std::sqrt(-dRho2));
                    cos_phi_.push_back(1.0);
                    sin_phi_.push_back(0.0);
                    cos_2phi_.push_back(1.0);
                    sin_2phi_.push_back(0.0);
                    break;

                default:
                    double dS = l - S_.back();

                    P_.push_back(0.0);
                    S_.push_back(dRho2 / dS - l);
                    cos_phi_.push_back(1.0);
                    sin_phi_.push_back(0.0);
                    cos_2phi_.push_back(1.0);
                    sin_2phi_.push_back(0.0);

                    break;
            }
        }
    }
}

double RadialWaveFunctions::channelPenetrationAndShift(unsigned int l, double rho2) const {
    double P = 0.0, S = 0.0;
    double rho = std::sqrt(std::abs(rho2));

    if (rho2 >= 0.0) {
        // Positive energy channel
        if (l == 0) {
            P = rho;
        } else {
            double P_prev = rho, S_prev = 0.0;
            for (unsigned int i = 1; i <= l; ++i) {
                double denom = std::pow(i - S_prev, 2) + std::pow(P_prev, 2);
                P = (rho2 * P_prev) / denom;
                S = (rho2 * (i - S_prev)) / denom - i;
                P_prev = P;
                S_prev = S;
            }
        }
    } else {
        // Negative energy channel
        if (l == 0) {
            S = -std::sqrt(-rho2);
        } else {
            double S_prev = -std::sqrt(-rho2);
            for (unsigned int i = 1; i <= l; ++i) {
                S = (rho2 / (i - S_prev)) - i;
                S_prev = S;
            }
        }
    }

    return P;
}

void RadialWaveFunctions::computePenetration(unsigned int l, double dRho2) const {
    double dP = P_.back();
    double dS = l - S_.back();
    double dC = dRho2 / (std::pow(dS, 2) + std::pow(dP, 2));
    P_.push_back(dC * dP);
}

void RadialWaveFunctions::computeShift(unsigned int l, double dRho2) const {
    double dP = P_.back();
    double dS = l - S_.back();
    double dC = dRho2 / (std::pow(dS, 2) + std::pow(dP, 2));
    S_.push_back(dC * dS - l);
}

void RadialWaveFunctions::computePhase(unsigned int l) const {
    double dP = P_.back();
    double dS = l - S_.back();
    double dY = dP / dS;
    double dN = 1.0 / (1.0 + std::pow(dY, 2));
    double dZ = std::sqrt(dN);
    double dCos = cos_phi_.back();
    double dSin = sin_phi_.back();
    double dCos_2 = cos_2phi_.back();
    double dSin_2 = sin_2phi_.back();

    cos_phi_.push_back((dCos + dY * dSin) * dZ);
    sin_phi_.push_back((dSin - dY * dCos) * dZ);
    cos_2phi_.push_back((1.0 - std::pow(dY, 2)) * dCos_2 + 2 * dY * dSin_2);
    sin_2phi_.push_back((1.0 - std::pow(dY, 2)) * dSin_2 - 2 * dY * dCos_2);
}