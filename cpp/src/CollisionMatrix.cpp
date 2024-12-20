#include "CollisionMatrix.h"
#include <complex>

CollisionMatrix::CollisionMatrix(const XMatrix& X, const SpinGroup& spinGroup) {
    const auto& channels = spinGroup.channels();
    size_t numChannels = channels.size();

    U_ = Eigen::MatrixXcd::Zero(numChannels, numChannels);

    const Eigen::MatrixXcd& X_matrix = X.getMatrix();

    for (size_t c = 0; c < numChannels; ++c) {
        double phi_c = channels[c].getPhaseShift();
        std::complex<double> phase_c = std::exp(-std::complex<double>(0, 1) * phi_c);

        for (size_t cp = 0; cp < numChannels; ++cp) {
            double phi_cp = channels[cp].getPhaseShift();
            std::complex<double> phase_cp = std::exp(std::complex<double>(0, 1) * phi_cp);

            if (c == cp) {
                U_(c, cp) = phase_c * (1.0 - 2.0 * std::complex<double>(0, 1) * X_matrix(c, cp)) * phase_cp;
            } else {
                U_(c, cp) = -2.0 * phase_c * std::complex<double>(0, 1) * X_matrix(c, cp) * phase_cp;
            }
        }
    }
}
