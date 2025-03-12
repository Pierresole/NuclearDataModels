#ifndef SPIN_GROUP_H
#define SPIN_GROUP_H

#include "Channel.h"
#include "ParticlePair.h"
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <iostream>

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
    // Track which channels are eliminated (true = eliminated, false = retained)
    std::vector<bool> eliminatedChannels_;

    // Method to return entrance ParticlePair
    ParticlePair entranceParticlePair() const {
        for (const auto& channel : channels_) {
            if (channel.getParticlePair().MT() == 2) {
                return channel.getParticlePair();
            }
        }
        throw std::runtime_error("No entrance channel (MT=2) found in spin group");
    };
    
    // Method to verify spin group attributes
    void verifyAttributes() const {
        // if (PJ_ != 1 && PJ_ != -1) {
        //     throw std::invalid_argument("PJ must be 1 or -1");
        // }
        
        // Verify each channel is valid
        for (const auto& channel : channels_) {
            if (!channel.isValid()) {
                throw std::invalid_argument("Invalid channel in spin group");
            }
        }
        
        // Verify each resonance has the correct number of gamma values
        for (const auto& res : resonances_) {
            if (res.getGamma().size() != channels_.size() && !channels_.empty()) {
                std::cout << "Number of gamma values in a resonance must match number of channels " << channels_.size() << std::endl;
                throw std::invalid_argument("Number of gamma values in a resonance must match number of channels");
            }
        }
    }

public:
    SpinGroup(double J, int PJ)
        : J_(J), PJ_(PJ) {
        verifyAttributes();
    }

    SpinGroup(double J, int PJ, const std::vector<Channel>& channels, const std::vector<Resonance>& resonances)
        : J_(J), PJ_(PJ), channels_(channels), resonances_(resonances) {
        // Initialize all channels as retained by default
        eliminatedChannels_.resize(channels_.size(), false);
        verifyAttributes();
    }

    void addChannel(const Channel& channel, bool isEliminated = false) { 
        channels_.push_back(channel); 
        // Add the elimination status for this channel
        eliminatedChannels_.push_back(isEliminated);
        // Verify the channel is valid
        if (!channel.isValid()) {
            channels_.pop_back();
            eliminatedChannels_.pop_back();
            throw std::invalid_argument("Invalid channel added to spin group");
        }
    }

    // Mark a channel as eliminated or retained
    void setChannelEliminated(size_t channelIndex, bool isEliminated = true) {
        if (channelIndex >= channels_.size()) {
            throw std::out_of_range("Channel index out of range");
        }
        eliminatedChannels_[channelIndex] = isEliminated;
    }

    // Check if a channel is eliminated
    bool isChannelEliminated(size_t channelIndex) const {
        if (channelIndex >= channels_.size()) {
            throw std::out_of_range("Channel index out of range");
        }
        return eliminatedChannels_[channelIndex];
    }

    // Get number of eliminated channels
    size_t getEliminatedChannelsCount() const {
        return std::count(eliminatedChannels_.begin(), eliminatedChannels_.end(), true);
    }

    // Get number of retained channels
    size_t getRetainedChannelsCount() const {
        return std::count(eliminatedChannels_.begin(), eliminatedChannels_.end(), false);
    }

    // Get indices of eliminated channels
    std::vector<size_t> getEliminatedChannelIndices() const {
        std::vector<size_t> indices;
        for (size_t i = 0; i < eliminatedChannels_.size(); ++i) {
            if (eliminatedChannels_[i]) {
                indices.push_back(i);
            }
        }
        return indices;
    }

    // Get indices of retained channels
    std::vector<size_t> getRetainedChannelIndices() const {
        std::vector<size_t> indices;
        for (size_t i = 0; i < eliminatedChannels_.size(); ++i) {
            if (!eliminatedChannels_[i]) {
                indices.push_back(i);
            }
        }
        return indices;
    }

    // Get eliminated channels
    std::vector<Channel> getEliminatedChannels() const {
        std::vector<Channel> result;
        for (size_t i = 0; i < channels_.size(); ++i) {
            if (eliminatedChannels_[i]) {
                result.push_back(channels_[i]);
            }
        }
        return result;
    }

    // Get retained channels
    std::vector<Channel> getRetainedChannels() const {
        std::vector<Channel> result;
        for (size_t i = 0; i < channels_.size(); ++i) {
            if (!eliminatedChannels_[i]) {
                result.push_back(channels_[i]);
            }
        }
        return result;
    }

    // Get specific eliminated channel
    const Channel& getEliminatedChannel(size_t eliminatedIndex) const {
        size_t count = 0;
        for (size_t i = 0; i < channels_.size(); ++i) {
            if (eliminatedChannels_[i]) {
                if (count == eliminatedIndex) {
                    return channels_[i];
                }
                count++;
            }
        }
        throw std::out_of_range("Eliminated channel index out of range");
    }

    // Get specific retained channel
    const Channel& getRetainedChannel(size_t retainedIndex) const {
        size_t count = 0;
        for (size_t i = 0; i < channels_.size(); ++i) {
            if (!eliminatedChannels_[i]) {
                if (count == retainedIndex) {
                    return channels_[i];
                }
                count++;
            }
        }
        throw std::out_of_range("Retained channel index out of range");
    }

    void addResonance(const Resonance& resonance) { 
        resonances_.push_back(resonance); 
        verifyAttributes();
    }

    // New overload that takes a resonance energy, widths, and conversion flag
    void addResonance(double energy, const std::vector<double>& widths, bool isReduced = true) {
        std::vector<double> reduced_widths;
        if (isReduced) {
            reduced_widths = widths; // Already reduced, use as is
        } else {
            reduced_widths.reserve(widths.size());
            // Convert from physical to reduced widths
            for (size_t i = 0; i < widths.size(); ++i) {
                reduced_widths.push_back(convertToReducedWidth(widths[i], energy, i));
            }
        }
        
        // Create and add the resonance with reduced widths
        Resonance resonance(energy, reduced_widths);
        resonances_.push_back(resonance);
        verifyAttributes();
    }

    const std::vector<Channel>& channels() const { return channels_; }
    const std::vector<Resonance>& getResonances() const { return resonances_; }
    double getJ() const { return J_; }
    int getPJ() const { return PJ_; }
    
    // Helper method to convert physical width to reduced width
    double convertToReducedWidth(double width, double energy, size_t channelIdx) const {
        if (width == 0.0) {
            return 0.0;
        }
        
        // Get the channel
        if (channelIdx >= channels_.size()) {
            throw std::out_of_range("Channel index out of range");
        }
        
        const Channel& channel = channels_[channelIdx];
        const ParticlePair& entrance = entranceParticlePair();
        
        // Check if this is a photon channel (similar to Python's check)
        // MT values 19 and 102 typically represent radiative capture
        double P = 1.0;
        if (channel.getParticlePair().MT() == 19 || channel.getParticlePair().MT() == 102) {
            // For photon channels, P = 1.0
            P = 1.0;
        } else {
            P = channel.computePenetration(energy, entrance);
            // P = channel.computeChannelQuantities(std::abs(energy), entrance).P;
        }

        
        if (P <= 0.0) {
            return 0.0;
        }
        
        // Convert to reduced width: γ = sign(Γ) * sqrt(|Γ|/(2*P))
        double reduced_width = (width >= 0 ? 1.0 : -1.0) * 
                               std::sqrt(std::abs(width) / (2.0 * P));
        
        return reduced_width;
    }

    // Helper method to check if this spin group is valid
    bool isValid() const {
        try {
            verifyAttributes();
            return true;
        } catch (const std::invalid_argument&) {
            return false;
        }
    }

    // Loop on channels of this spingroup
    double crossSection(double E, const ParticlePair& entrancePP) const {
        Eigen::MatrixXcd U = computeCollisionMatrix(E, entrancePP);
        
        double sigma = 0.0;
        int nChannels = channels_.size();
        for (int c = 0; c < nChannels; ++c) {
            double partialSigmaJ = (2*M_PI/(entrancePP.k2(E, entrancePP))) * (1. - U(c,c).real());
            sigma += partialSigmaJ;
        }
        return sigma;
    }

    Eigen::MatrixXcd computeCollisionMatrix(double E, const ParticlePair& entrancePP) const {
        int nChannels = channels_.size();
        Eigen::MatrixXcd Collision(nChannels,nChannels);
        
        std::vector<Channel::ChannelQuantities> channelQuantities;
        channelQuantities.reserve(nChannels);
        
        // First calculate all channel quantities and store them
        for (int c = 0; c < nChannels; ++c) {
            channelQuantities.push_back(channels_[c].computeChannelQuantities(E, entrancePP));
        }

        // Pass pre-computed quantities to computeXMatrix
        auto X = computeXMatrix(E, channelQuantities);

        for (int c = 0; c < nChannels; ++c) {
            const auto& cQuant = channelQuantities[c];
            double sqrtP_c = sqrt(cQuant.P);
            std::complex<double> phase_c = std::exp(-std::complex<double>(0, 1) * cQuant.phi);

            for (int cp = 0; cp < nChannels; ++cp) {
                const auto& cpQuant = channelQuantities[cp];
                double sqrtP_cp = sqrt(cpQuant.P);
                std::complex<double> phase_cp = std::exp( - std::complex<double>(0, 1) * cpQuant.phi);

                // if (c == cp) {
                    Collision(c, cp) = phase_c * (1.0 + 2.0 * std::complex<double>(0, 1) * X(c, cp)) * phase_cp;
                // } else {
                //     Collision(c, cp) = -2.0 * phase_c * std::complex<double>(0, 1) * X(c, cp) * phase_cp;
                // }

                // If elastic
                // U(c,cp) = exp(-2i*phi_c)*( (c==cp ? 1.0 : 0.0) + 2i*sqrtP_c*X(c,cp)*sqrtP_cp );
            }
        }

        return Collision;
    }

    // Modified to handle eliminated and retained channels
    Eigen::MatrixXcd computeXMatrix(double E, const std::vector<Channel::ChannelQuantities>& channelQuants) const {
        // Get total number of channels and retained channels
        int nTotalChannels = channels_.size();
        std::vector<size_t> retainedIndices = getRetainedChannelIndices();
        std::vector<size_t> eliminatedIndices = getEliminatedChannelIndices();
        int nRetained = retainedIndices.size();
        int nEliminated = eliminatedIndices.size();
        
        // Create reduced R matrix directly for retained channels
        Eigen::MatrixXcd Rcc = Eigen::MatrixXcd::Zero(nRetained, nRetained);
        
        // Populate R-matrix for retained channels with modified denominator
        for (int i = 0; i < nRetained; ++i) {
            for (int j = 0; j < nRetained; ++j) {
                size_t c = retainedIndices[i];
                size_t cp = retainedIndices[j];
                
                // Calculate R matrix element for this pair of retained channels
                for (const auto& resonance : resonances_) {
                    double gamma_c = resonance.getGamma()[c];
                    double gamma_cp = resonance.getGamma()[cp];
                    
                    // Base denominator
                    std::complex<double> denominator(resonance.getEnergy() - E, 0.0);
                    
                    // Add contribution from eliminated channels to denominator
                    if (nEliminated > 0) {
                        std::complex<double> D(0.0, 0.0);
                        
                        // Sum over all eliminated channels
                        for (int e = 0; e < nEliminated; ++e) {
                            size_t eliminatedIdx = eliminatedIndices[e];
                            double gamma_e = resonance.getGamma()[eliminatedIdx];
                            
                            // Get channel quantities for this eliminated channel
                            const auto& quant = channelQuants[eliminatedIdx];
                            std::complex<double> L_e = quant.S + std::complex<double>(0, quant.P);
                            // Accumulate contribution
                            D += std::pow(gamma_e, 2) * L_e;
                        }
                        
                        // Subtract the contribution from eliminated channels
                        denominator -= D;
                    }
                    
                    // Add this resonance's contribution to the R matrix element
                    Rcc(i, j) += gamma_c * gamma_cp / denominator;
                }
            }
        }
        
        // Now create the X matrix using only retained channels
        Eigen::MatrixXcd X(nRetained, nRetained);
        Eigen::MatrixXcd L(nRetained, nRetained), sqrtP(nRetained, nRetained);
        
        // Fill the L and sqrtP matrices for retained channels
        for (int i = 0; i < nRetained; ++i) {
            const auto& quant = channelQuants[retainedIndices[i]];
            L(i, i) = quant.S + std::complex<double>(0, quant.P);
            sqrtP(i, i) = std::sqrt(quant.P);
        }
        
        // Calculate X matrix for retained channels
        X = sqrtP * (Eigen::MatrixXcd::Identity(nRetained, nRetained) - Rcc * L).inverse() * Rcc * sqrtP;
        
        // std::cout << "========================" << std::endl;
        // std::cout << "R" << std::endl;
        // std::cout << Rcc << std::endl;
        // std::cout << "X" << std::endl;
        // std::cout << X << std::endl;
        // std::cout << "Id - R*L" << std::endl;
        // std::cout << (Eigen::MatrixXcd::Identity(nRetained, nRetained) - Rcc * L) << std::endl;
        // std::cout << "(Id - R*L)^-1 * R" << std::endl;
        // std::cout << (Eigen::MatrixXcd::Identity(nRetained, nRetained) - Rcc * L).inverse() * Rcc << std::endl;

        // Create and return full X matrix with zeros for eliminated channels
        Eigen::MatrixXcd fullX = Eigen::MatrixXcd::Zero(nTotalChannels, nTotalChannels);
        for (int i = 0; i < nRetained; ++i) {
            for (int j = 0; j < nRetained; ++j) {
                fullX(retainedIndices[i], retainedIndices[j]) = X(i, j);
            }
        }
        
        return fullX;
    }
    
    void FillRMatrix(Eigen::MatrixXcd& R, double E) const {
        const auto& resonances = this->getResonances();
        const auto& channels = this->channels();
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
                R(c, cp) = sum;
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
                std::complex<double> sum = 0.0;
                int iChCapt = 0;
                for (const auto& res : resonances) {
                    sum += res.getGamma()[c] * res.getGamma()[cp] / (res.getEnergy() - E - std::complex<double>(0, res.getGamma()[iChCapt]/2));
                }
                R(c, cp) = sum;
            }
        }
    }
};

#endif // SPIN_GROUP_H
