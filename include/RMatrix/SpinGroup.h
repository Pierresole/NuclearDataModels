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

    double computeCrossSection(double E) const {
        Eigen::MatrixXd R = computeR(E);
        Eigen::MatrixXcd X = computeX(R, E);
        CollisionMatrix U(X, channels_, E);

        double sigma_group = 0.0;
        for (size_t c = 0; c < channels_.size(); ++c) {
            sigma += computeChannelCrossSection(U, c, E);
        }
        return sigma;
    }

    
    double crossSection(double E, const ParticlePair& particlePair) const {
        Eigen::MatrixXcd U = computeCollisionMatrix(E);
        double k = particlePair.waveNumber(E);
        double sigma = 0.0;
        int nCh = channels_.size();
        for (int c = 0; c < nCh; ++c) {
            // typically elastic cross section: c=c'
            sigma += (M_PI/(k*k)) * std::norm(U(c,c) - 1.);
        }
        return sigma;
    }
    

    Eigen::MatrixXd FillRMatrix(const SpinGroup& spinGroup, Eigen::MatrixXcd& R, double E) {
        const auto& resonances = spinGroup.getResonances();
        const auto& channels = spinGroup.channels();
        size_t numChannels = channels.size();
    
        // R_ = Eigen::MatrixXd::Zero(numChannels, numChannels);
    
        for (size_t c = 0; c < numChannels; ++c) {
            for (size_t cp = 0; cp < numChannels; ++cp) {
                double sum = 0.0;
                for (const auto& resonance : resonances) {
                    double gamma_c = resonance.getGamma()[c];
                    double gamma_cp = resonance.getGamma()[cp];
                    sum += gamma_c * gamma_cp / (resonance.getEnergy() - E);
                }
                R_(c, cp) = sum;
            }
        }
    }

    Eigen::MatrixXd FillReichMooreRMatrix(const SpinGroup& spinGroup, Eigen::MatrixXcd& R, double E) {
        const auto& resonances = spinGroup.getResonances();
        const auto& channels = spinGroup.channels();
        size_t numChannels = channels.size();
    
        // R_ = Eigen::MatrixXd::Zero(numChannels, numChannels);
    
        for (size_t c = 0; c < numChannels; ++c) {
            for (size_t cp = 0; cp < numChannels; ++cp) {
                double sum = 0.0;
                for (const auto& resonance : resonances) {
                    sum += res.gamma[c] * res.gamma[cp] / (res.energy - E - std::complex<double>(0, res.gammaGamma/2));
                }
                R_(c, cp) = sum;
            }
        }
    }

    Eigen::MatrixXcd computeXMatrix(double E) const {
        int nCh = channels_.size();
        Eigen::MatrixXcd R(nCh, nCh), L(nCh, nCh), sqrtP(nCh, nCh);
    
        // R-matrix calculation (Reich-Moore)
        FillRMatrix(*this, R, E);
    
        for (int c = 0; c < nCh; ++c) {
            L(c, c) = channels_[c].shiftFactor(E) + std::complex<double>(0, channels_[c].penetrability(E));
            sqrtP(c, c) = sqrt(channels_[c].penetrability(E));
        }
    
        // Check if E matches a resonance energy
        bool isResonanceEnergy = false;
        double epsilon = 1e-6; // tolerance for numerical equality
        for (const auto& res : resonances_) {
            if (std::abs(E - res.energy) < epsilon) {
                isResonanceEnergy = true;
                break;
            }
        }
    
        Eigen::MatrixXcd X(nCh, nCh);
    
        if (isResonanceEnergy) {
            // Use ST decomposition
            Eigen::MatrixXcd S = Eigen::MatrixXcd::Identity(nCh, nCh);
            Eigen::MatrixXcd T = R;
            double closestEnergy = E; // since E matches resonance
    
            for (int c = 0; c < nCh; ++c) {
                S(c, c) *= (closestEnergy - E); // Here, it's zero, so directly T is adjusted
                T.row(c) *= (closestEnergy - E);
            }
    
            X = sqrtP * (S - T * L).inverse() * T * sqrtP;
        } else {
            // Regular calculation
            X = sqrtP * (Eigen::MatrixXcd::Identity(nCh,nCh) - R*L).inverse() * R * sqrtP;
        }
    
        return X;
    }
    
    Eigen::MatrixXcd computeCollisionMatrixU(double E) const {
        auto X = computeXMatrix(E);
        int nCh = channels_.size();
        Eigen::MatrixXcd U(nCh,nCh);
        for (int c = 0; c < nCh; ++c) {
            double phi = channels_[c].phaseShift(E);
            double sqrtP_c = sqrt(channels_[c].penetrability(E));
            for (int cp = 0; cp < nCh; ++cp) {
                double sqrtP_cp = sqrt(channels_[cp].penetrability(E));
                U(c,cp) = exp(-2i*phi)*( (c==cp ? 1.0 : 0.0) + 2i*sqrtP_c*X(c,cp)*sqrtP_cp );
            }
        }
        return U;
    }


#endif // SPIN_GROUP_H
