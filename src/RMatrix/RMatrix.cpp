#include "RMatrix.h"

RMatrix::RMatrix(const SpinGroup& spinGroup, double E) {
    const auto& resonances = spinGroup.getResonances();
    const auto& channels = spinGroup.channels();
    size_t numChannels = channels.size();

    R_ = Eigen::MatrixXd::Zero(numChannels, numChannels);

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
