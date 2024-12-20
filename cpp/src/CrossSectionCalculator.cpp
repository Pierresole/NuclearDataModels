#include "CrossSectionCalculator.h"
#include <iostream>
#include <cmath>

CrossSectionCalculator::CrossSectionCalculator(const std::vector<SpinGroup>& spinGroups, double A)
    : spinGroups_(spinGroups), A_(A)
{
}

bool CrossSectionCalculator::channelMatches(const Channel& channel, const std::string& reactionChannel) const {
    // Simple comparison based on particle pair identifiers
    return channel.getParticlePair().getReactionID() == reactionChannel;
}

std::vector<double> CrossSectionCalculator::computeCrossSections(const std::vector<double>& energies, const std::string& reactionChannel) {
    std::vector<double> crossSections; // To store total cross sections for each energy

    for (double E : energies) {
        double totalCrossSection = 0.0;
        std::cout << "E=" << E << std::endl;
        for (const auto& spingroup : spinGroups_) {
            std::cout << "Spingroup=" << spingroup.getJ() << "^" << spingroup.getPJ() << "  (/" << spinGroups_.size() << ")" << std::endl;
            // Filter channels based on reactionChannel
            std::vector<Channel> selectedChannels;
            for (const auto& channel : spingroup.channels()) {
                if (reactionChannel == "(n,tot)" || channelMatches(channel, reactionChannel)) {
                    selectedChannels.push_back(channel);
                    std::cout << "Channel " << channel.getParticlePair().getReactionID() << " selected" << std::endl;
                }
            }

            if (selectedChannels.empty()) {
                continue; // No matching channels in this spin group
            }

            // Compute channel-specific quantities for selected channels
            for (auto& channel : selectedChannels) {
                //channel.computeChannelQuantities(E);
            }

            // Build R-matrix
            RMatrix R(spingroup, E);

            // Build P and S matrices for selected channels
            size_t numChannels = selectedChannels.size();
            Eigen::MatrixXd P = Eigen::MatrixXd::Zero(numChannels, numChannels);
            Eigen::MatrixXd S = Eigen::MatrixXd::Zero(numChannels, numChannels);

            std::cout << "Computing Penetration and Shift..." << std::endl;
            for (size_t c = 0; c < numChannels; ++c) {
                P(c, c) = selectedChannels[c].getPenetrationFactor();
                S(c, c) = selectedChannels[c].getShiftFactor();
            }

            std::cout << "P = \n" << P << std::endl;
            std::cout << "S = \n" << S << std::endl;

            // Compute L matrix
            std::cout << "Computing L=S+iP..." << std::endl;
            Eigen::MatrixXcd L = S.cast<std::complex<double>>() + std::complex<double>(0, 1) * P.cast<std::complex<double>>();
            std::cout << "L = \n" << L << std::endl;

            std::cout << "R size " << R.getMatrix().rows() << " " << R.getMatrix().cols() << std::endl;
            std::cout << "L size " << L.rows() << " " << L.cols() << std::endl;

            // Compute Level Matrix A
            std::cout << "Computing A=(I-RL)^-1..." << std::endl;
            LevelMatrix A(R.getMatrix().cast<std::complex<double>>(), L);

            // LevelMatrix A(R.getMatrix(), L);

            // // Compute X Matrix
            // std::cout << "Computing X=PcPc'SumA..." << std::endl;
            // XMatrix X(A, spingroup, E);

            // // Compute Collision Matrix U
            // std::cout << "Computing U=eiPhi(1+2iX)..." << std::endl;
            // CollisionMatrix U(X, spingroup);

            // // Compute cross sections for selected channels
            // const Eigen::MatrixXcd& U_matrix = U.getMatrix();

            // for (size_t c = 0; c < numChannels; ++c) {
            //     double k_c = selectedChannels[c].getWaveNumber();

            //     if (reactionChannel == "(n,n)" || reactionChannel == "(n,tot)") {
            //         // Elastic Cross Section
            //         double sigma_elastic = (M_PI / (k_c * k_c)) * std::norm(1.0 - U_matrix(c, c));

            //         // Accumulate the total cross section
            //         totalCrossSection += sigma_elastic;
            //     }

            //     for (size_t cp = 0; cp < numChannels; ++cp) {
            //         if (c != cp) {
            //             if (channelMatches(selectedChannels[cp], reactionChannel) || reactionChannel == "(n,tot)") {
            //                 // Reaction Cross Section
            //                 double sigma_reaction = (M_PI / (k_c * k_c)) * std::norm(U_matrix(c, cp));

            //                 // Accumulate the total cross section
            //                 totalCrossSection += sigma_reaction;
            //             }
            //         }
            //     }
            // }
        }

        // Store the total cross section for this energy
        crossSections.push_back(totalCrossSection);
    }

    return crossSections; // Return the vector of cross sections
}