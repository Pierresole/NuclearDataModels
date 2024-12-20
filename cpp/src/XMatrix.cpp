#include "XMatrix.h"
#include <cmath>

XMatrix::XMatrix(const LevelMatrix& A, const SpinGroup& spinGroup, double E) {
    const auto& resonances = spinGroup.getResonances();
    const auto& channels = spinGroup.channels();
    size_t numChannels = channels.size();
    size_t numResonances = resonances.size();

    Eigen::MatrixXd Gamma = Eigen::MatrixXd::Zero(numResonances, numChannels);

    // Compute partial widths
    for (size_t n = 0; n < numResonances; ++n) {
        const auto& gamma = resonances[n].getGamma();
        for (size_t c = 0; c < numChannels; ++c) {
            double gamma_nc = gamma[c];
            double P_c = channels[c].getPenetrationFactor();
            Gamma(n, c) = 2.0 * P_c * gamma_nc * gamma_nc;
        }
    }

    // Construct X-matrix
    X_ = Eigen::MatrixXcd::Zero(numChannels, numChannels);

    const Eigen::MatrixXcd& A_matrix = A.getMatrix();

    for (size_t c = 0; c < numChannels; ++c) {
        for (size_t cp = 0; cp < numChannels; ++cp) {
            std::complex<double> sum = 0.0;
            for (size_t lambda = 0; lambda < numResonances; ++lambda) {
                for (size_t mu = 0; mu < numResonances; ++mu) {
                    double sqrt_Gamma_lambda_c = std::sqrt(Gamma(lambda, c));
                    double sqrt_Gamma_mu_cp = std::sqrt(Gamma(mu, cp));
                    sum += sqrt_Gamma_lambda_c * A_matrix(lambda, mu) * sqrt_Gamma_mu_cp;
                }
            }
            X_(c, cp) = 0.5 * sum;
        }
    }
}
