// SAlphaBeta.cpp

#include "ThermalScattering.h"

SAlphaBeta::SAlphaBeta(
    double temperature,
    double w_s,
    const std::vector<double>& alpha_grid,
    const std::vector<double>& beta_grid,
    int n_max,
    const std::vector<double>& energy_grid,
    const std::vector<double>& rho_energy,
    const std::vector<std::tuple<double, double, std::string, double>>& peaks)
    : temperature_(temperature), w_s_(w_s), alpha_grid_(alpha_grid), beta_grid_(beta_grid), n_max_(n_max), energy_grid_(energy_grid), rho_energy_(rho_energy)
{
    // Initialize member variables
    lambda_ = 0.0;
    T_s_bar_ = 0.0;
    alpha_max_.resize(beta_grid_.size());
    S_alpha_beta_.resize(beta_grid_.size(), std::vector<double>(alpha_grid_.size(), 0.0));

    // Initialize P_beta_ and T1_beta_
    P_beta_.resize(beta_grid_.size(), 0.0);
    T1_beta_.resize(beta_grid_.size(), 0.0);

    // Step 1: Read and Convert the Phonon Spectrum
    read_and_convert_spectrum();

    // Compute S(alpha, beta)
    compute_S_alpha_beta();
}

// Helper function for linear interpolation
double interpolate(double x, const std::vector<double>& x_vals, const std::vector<double>& y_vals)
{
    if (x <= x_vals.front()) {
        return y_vals.front();
    }
    if (x >= x_vals.back()) {
        return y_vals.back();
    }

    auto it = std::lower_bound(x_vals.begin(), x_vals.end(), x);
    size_t idx = std::distance(x_vals.begin(), it);

    // Ensure idx is at least 1 to avoid accessing negative indices
    if (idx == 0) idx = 1;

    double x0 = x_vals[idx - 1];
    double x1 = x_vals[idx];
    double y0 = y_vals[idx - 1];
    double y1 = y_vals[idx];

    double t = (x - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
}


void SAlphaBeta::read_and_convert_spectrum()
{
    // Step 1: Define beta_spectrum_grid_
    size_t num_points = rho_energy_.size();
    beta_spectrum_grid_.resize(num_points);

    // Define beta_spectrum_grid_ as uniformly spaced beta values
    double beta_min = energy_grid_.front() / (k_B * temperature_);
    double beta_max = energy_grid_.back() / (k_B * temperature_);

    for (size_t i = 0; i < num_points; ++i) {
        beta_spectrum_grid_[i] = beta_min + i * (beta_max - beta_min) / (num_points - 1);
    }

    // Step 2: Interpolate rho_energy_ at energies corresponding to beta values
    rho_beta_.resize(num_points);

    for (size_t i = 0; i < num_points; ++i) {
        double beta = beta_spectrum_grid_[i];
        double e = beta * k_B * temperature_;
        // Interpolate rho_energy_ at energy e
        rho_beta_[i] = interpolate(e, energy_grid_, rho_energy_);
    }

    // Step 3: Adjust rho_beta_ at beta = 0 to enforce beta squared behavior
    // if (beta_spectrum_grid_.size() > 1) {
    //     rho_beta_[0] = rho_beta_[1] / std::pow(beta_spectrum_grid_[1], 2);
    // }

    // Step 4: Normalize rho(beta) over beta_spectrum_grid_
    double integral_rho_beta = 0.0;
    for (size_t i = 0; i < beta_spectrum_grid_.size() - 1; ++i) {
        double delta_beta = beta_spectrum_grid_[i + 1] - beta_spectrum_grid_[i];
        double avg_rho = (rho_beta_[i] + rho_beta_[i + 1]) / 2.0;
        integral_rho_beta += avg_rho * delta_beta;
    }

    if (integral_rho_beta > 0.0) {
        double normalization_factor = w_s_ / integral_rho_beta;
        for (size_t i = 0; i < rho_beta_.size(); ++i) {
            rho_beta_[i] *= normalization_factor;
        }
    } else {
        std::cerr << "Warning: Integral of rho(beta) is zero or negative." << std::endl;
    }
}

// Linear interpolation of rho(beta)
std::function<double(double)> interpolate_rho_beta(const std::vector<double>& beta_spectrum_grid_,  const std::vector<double>& rho_beta_) 
{
    return [beta_spectrum_grid_, rho_beta_](double beta) {
        // Ensure beta is within bounds
        if (beta <= beta_spectrum_grid_.front()) {
            return rho_beta_.front();
        }
        if (beta >= beta_spectrum_grid_.back()) {
            return rho_beta_.back();
        }

        // Find the interval in beta_spectrum_grid_ that contains beta
        auto upper = std::upper_bound(beta_spectrum_grid_.begin(), beta_spectrum_grid_.end(), beta);
        size_t i = std::distance(beta_spectrum_grid_.begin(), upper) - 1;

        // Linear interpolation
        double beta0 = beta_spectrum_grid_[i];
        double beta1 = beta_spectrum_grid_[i + 1];
        double rho0 = rho_beta_[i];
        double rho1 = rho_beta_[i + 1];

        return rho0 + (beta - beta0) * (rho1 - rho0) / (beta1 - beta0);
    };
}

double SAlphaBeta::compute_P_beta(double beta, const std::function<double(double)>& rho_beta_func) {
    if (beta == 0) {
        return rho_beta_func(beta_spectrum_grid_[1]) / std::pow(beta_spectrum_grid_[1], 2);
    }

    // Handle linear interpolation for beta between 0 and beta_spectrum_grid_[1]
    if (beta > 0 && beta < beta_spectrum_grid_[1]) {
        double P_beta_0 = rho_beta_func(beta_spectrum_grid_[1]) / std::pow(beta_spectrum_grid_[1], 2);  // P(0)
        double P_beta_1 = rho_beta_func(beta_spectrum_grid_[1]) / (2 * beta_spectrum_grid_[1] * std::sinh(beta_spectrum_grid_[1] / 2.0));  // P(beta_spectrum_grid_[1])

        // Linear interpolation
        return P_beta_0 + (P_beta_1 - P_beta_0) * (beta / beta_spectrum_grid_[1]);
    }

    // Default case for beta >= beta_spectrum_grid_[1]
    return rho_beta_func(beta) / (2 * beta * std::sinh(beta / 2.0));
}

double SAlphaBeta::compute_lambda(const std::function<double(double)>& P_beta_func) {
    //! Formula implemented :
    //!   \int_{-\infty}^{+\infty} P(B) * exp(-B / 2) dB
    //! = \int_{   0   }^{+\infty} P(B) * exp(-B / 2) dB + \int_{-\infty}^{   0   } P(B)   * exp(-B / 2) dB
    //! = \int_{   0   }^{+\infty} P(B) * exp(-B / 2) dB + \int_{   0   }^{+\infty} P(-B') * exp(-B' / 2) dB' with B' = - B
    //! = \int_{   0   }^{+\infty} P(B) * [ exp(-B / 2) + exp(B / 2)] dB    assuming P(-B)=P(B) 
    //! = \int_{   0   }^{+\infty} P(B) * 2 * cosh(B / 2) 
    
    return boost::math::quadrature::trapezoidal([P_beta_func](double B) { return P_beta_func(B) * 2 * std::cosh(B / 2.0); }, 0.0, beta_spectrum_grid_.back(), 1e-8);
}

void SAlphaBeta::compute_effective_temperature(const std::function<double(double)>& P_beta_func)
{
    //! Formula implemented :
    //!   T / 2w_s * \int_{-\infty}^{+\infty} B^2 * P(B) * exp(-B) dB
    //!   (trick [0->+\infty] + [-\infty->0])
    //!    =  T / 2w_s * \int_{0}^{+\infty} B^2 * P(B) * 2 * cosh(B) dB
    //T_s_bar_ = (temperature_ / w_s_) * boost::math::quadrature::trapezoidal([P_beta_func](double B) { return std::pow(B,2) * P_beta_func(B) * std::cosh(B); }, 0.0, beta_spectrum_grid_.back(), 1e-8);
    // NJOY : cosh(beta/2)
    T_s_bar_ = (temperature_ / w_s_) * boost::math::quadrature::trapezoidal([P_beta_func](double B) { return std::pow(B,2) * P_beta_func(B) * std::cosh(0.5 * B); }, 0.0, beta_spectrum_grid_.back(), 1e-8);
}

double SAlphaBeta::compute_T1(double beta, const std::function<double(double)>& P_beta_func) {
    return P_beta_func(beta) * exp( -beta / 2.0) / lambda_;
}

double compute_Tn(int n,  double beta,  const std::function<double(double)>& T1_func,  const std::function<double(int, double)>& Tn_func) 
{
    if (n == 1) return T1_func(beta);

    // Use Boost quadrature for the integral
    boost::math::quadrature::tanh_sinh<double> integrator;
    double result = integrator.integrate(
        [T1_func, Tn_func, n, beta](double beta_prime) {return T1_func(beta - beta_prime) * Tn_func(n - 1, beta_prime);}, 
        0.0, beta);

    return result;
}

double compute_SAB(double alpha, double beta, int n_max, double lambda, const std::function<double(int, double)>& Tn_func) 
{
    // Exponential term
    double exp_term = std::exp(-alpha * lambda);
    
    // Sum term
    double sum_term = 0.0;
    for (int n = 1; n <= n_max; ++n) {
        double alpha_lambda_pow_n = std::pow(alpha * lambda, n);
        double factorial_n = std::tgamma(n + 1);  // n!
        double Tn_beta = Tn_func(n, beta);
        sum_term += (alpha_lambda_pow_n / factorial_n) * Tn_beta;
    }

    return exp_term * sum_term;
}


std::vector<double> SAlphaBeta::convolve(const std::vector<double>& f, const std::vector<double>& g)
{
    // Perform convolution using numerical integration over the beta grid
    size_t N = beta_grid_.size();
    std::vector<double> result(N, 0.0);
    double beta_min_ = beta_grid_.front(); 
    double beta_max_ = beta_grid_.back(); 

    for (size_t i = 0; i < N; ++i) {
        double beta = beta_grid_[i];

        // Integrate over beta'
        double integral = boost::math::quadrature::trapezoidal(
            [&](double beta_prime) {
                double beta_minus_beta_prime = beta - beta_prime;

                // Interpolate f and g at beta_prime and beta - beta_prime
                double f_val = 0.0, g_val = 0.0;

                if (beta_prime >= beta_min_ && beta_prime <= beta_max_) {
                    auto it_f = lower_bound(beta_grid_.begin(), beta_grid_.end(), beta_prime);
                    size_t idx_f = distance(beta_grid_.begin(), it_f);
                    if (idx_f >= N) idx_f = N - 1;
                    f_val = f[idx_f];
                }

                if (beta_minus_beta_prime >= beta_min_ && beta_minus_beta_prime <= beta_max_) {
                    auto it_g = lower_bound(beta_grid_.begin(), beta_grid_.end(), beta_minus_beta_prime);
                    size_t idx_g = distance(beta_grid_.begin(), it_g);
                    if (idx_g >= N) idx_g = N - 1;
                    g_val = g[idx_g];
                }

                return f_val * g_val;
            },
            beta_min_, beta_max_, 1e-6);

        result[i] = integral;
    }

    return result;
}

double SAlphaBeta::interpolate_ln_T(const std::vector<double>& ln_T, double beta)
{
    // Check if beta is outside the grid range
    if (beta <= beta_grid_.front()) {
        return ln_T.front();
    }
    if (beta >= beta_grid_.back()) {
        return ln_T.back();
    }

    // Find the interval containing beta
    auto it = std::lower_bound(beta_grid_.begin(), beta_grid_.end(), beta);
    size_t idx = std::distance(beta_grid_.begin(), it);

    // Linear interpolation in log space
    double beta1 = beta_grid_[idx - 1];
    double beta2 = beta_grid_[idx];
    double ln_T1 = ln_T[idx - 1];
    double ln_T2 = ln_T[idx];

    double slope = (ln_T2 - ln_T1) / (beta2 - beta1);
    double ln_T_val = ln_T1 + slope * (beta - beta1);

    return ln_T_val;
}


std::vector<double> SAlphaBeta::compute_ln_Tn(const std::vector<double>& ln_T1, const std::vector<double>& ln_Tn_minus1)
{
    size_t N = beta_grid_.size();
    std::vector<double> ln_Tn(N, -std::numeric_limits<double>::infinity()); // Initialize with log(0)

    // For each beta
    for (size_t i = 0; i < N; ++i) {
        double beta = beta_grid_[i];

        // Define the integrand for convolution
        auto integrand = [&](double beta_prime) -> double {
            double beta_minus_beta_prime = beta - beta_prime;

            // Ensure beta_prime and beta_minus_beta_prime are within the grid
            if (beta_prime < beta_grid_.front() || beta_prime > beta_grid_.back() ||
                beta_minus_beta_prime < beta_grid_.front() || beta_minus_beta_prime > beta_grid_.back()) {
                return 0.0; // Return 0 for out-of-range values
            }

            // Interpolate ln_T1 and ln_Tn_minus1
            double ln_T1_val = interpolate_ln_T(ln_T1, beta_prime);
            double ln_Tn_minus1_val = interpolate_ln_T(ln_Tn_minus1, beta_minus_beta_prime);

            // Sum the logarithms to multiply the values (since ln(a*b) = ln(a) + ln(b))
            double ln_product = ln_T1_val + ln_Tn_minus1_val;

            // Exponentiate to get the normal-space value for integration
            return std::exp(ln_product);
        };

        // Perform the integration over beta_prime
        double integral = boost::math::quadrature::trapezoidal(
            integrand, beta_grid_.front(), beta_grid_.back(), 1e-6);

        // Compute ln(Tn(beta))
        if (integral > 0.0) {
            ln_Tn[i] = std::log(integral);
        } else {
            ln_Tn[i] = -std::numeric_limits<double>::infinity(); // log(0)
        }
    }

    return ln_Tn;
}

double log_sum_exp(double x, double y)
{
    if (x == -std::numeric_limits<double>::infinity()) return y;
    if (y == -std::numeric_limits<double>::infinity()) return x;
    double max_val = std::max(x, y);
    return max_val + std::log1p(std::exp(-std::abs(x - y)));
}

/*
void SAlphaBeta::compute_S_alpha_beta()
{
    size_t n_alpha = alpha_grid_.size();
    size_t n_beta = beta_grid_.size();

    // Create an interpolated rho(beta) function
    auto rho_beta_func = interpolate_rho_beta(beta_spectrum_grid_, rho_beta_);
    
    // Compute lambda
    auto P = [this, rho_beta_func](double beta) { return compute_P_beta(beta, rho_beta_func); };
    lambda_ = compute_lambda(P);
    compute_effective_temperature(P);
    auto T1_func = [this, P](double beta) { return compute_T1(beta, P); };

    output_debug_information();

    std::vector<double> ln_T1;
    for(double beta : beta_grid_) {ln_T1.push_back( std::log(T1_func(beta)) );}


    // For each alpha
    for (size_t i_alpha = 0; i_alpha < n_alpha; ++i_alpha) {
        double alpha = alpha_grid_[i_alpha];

        // Initialize ln_Tn_list to store ln Tn(beta)
        std::vector<std::vector<double>> ln_Tn_list;
        ln_Tn_list.push_back(ln_T1); // T1 corresponds to n = 1

        // Compute ln of the exponential factor
        double ln_exp_factor = -alpha * lambda_;

        // Initialize ln S(alpha, beta) with log(0)
        std::vector<double> ln_S_beta(n_beta, -std::numeric_limits<double>::infinity());

        // Initialize ln_coef for n = 1
        double ln_coef = std::log(alpha * lambda_) - std::log(1); // ln((alpha * lambda)^1 / 1!)

        // Summation over n
        for (int n = 1; n <= n_max_; ++n) {
            if (n > 1) {
                ln_coef += std::log(alpha * lambda_) - std::log(n);
            }

            // Compute ln Tn(beta)
            std::vector<double> ln_Tn;
            if (n == 1) {
                ln_Tn = ln_T1;
            } else {
                ln_Tn = compute_ln_Tn(ln_T1, ln_Tn_list[n - 2]);
                ln_Tn_list.push_back(ln_Tn);
            }

            // Add the term to ln S_beta
            for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
                double ln_term = ln_coef + ln_Tn[i_beta];

                // Update ln S_beta using log-sum-exp
                ln_S_beta[i_beta] = log_sum_exp(ln_S_beta[i_beta], ln_term);
            }

            // **Convergence Criteria**

            // Check if the maximum ln_term is significantly smaller than ln S_beta
            double max_ln_term = *std::max_element(ln_Tn.begin(), ln_Tn.end()) + ln_coef;
            double max_ln_S_beta = *std::max_element(ln_S_beta.begin(), ln_S_beta.end());

            // Difference between the new term and current sum
            double diff = max_ln_term - max_ln_S_beta;

            // Set a threshold (e.g., contributions less than 1e-8)
            if (diff < std::log(1e-4)) {
                break; // Convergence achieved
            }
        }

        // Apply exponential factor to ln S_beta
        for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
            ln_S_beta[i_beta] += ln_exp_factor;
        }

        // Convert ln S_beta to S_beta
        std::vector<double> S_beta(n_beta, 0.0);
        for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
            S_beta[i_beta] = std::exp(ln_S_beta[i_beta]);
        }

        // Store S_beta in S_alpha_beta_
        S_alpha_beta_[i_alpha] = S_beta;
    }
}
*/

/*
void SAlphaBeta::compute_S_alpha_beta()
{
    size_t n_alpha = alpha_grid_.size();
    size_t n_beta = beta_grid_.size();

    // Create an interpolated rho(beta) function
    auto rho_beta_func = interpolate_rho_beta(beta_spectrum_grid_, rho_beta_);
    
    // Compute lambda
    auto P = [this, rho_beta_func](double beta) { return compute_P_beta(beta, rho_beta_func); };
    lambda_ = compute_lambda(P);
    compute_effective_temperature(P);
    auto T1_func = [this, P](double beta) { return compute_T1(beta, P); };

    output_debug_information();

    std::vector<double> T1;
    for(double beta : beta_grid_) {
        T1.push_back(T1_func(beta));
    }

    // For each alpha
    for (size_t i_alpha = 0; i_alpha < n_alpha; ++i_alpha) {
        double alpha = alpha_grid_[i_alpha];

        // Initialize Tn(beta)
        std::vector<std::vector<double>> Tn_list;
        Tn_list.push_back(T1);

        // Compute exponential factor
        double exp_factor = exp(-alpha * lambda_);

        // Initialize S(alpha, beta) for this alpha
        std::vector<double> S_beta(n_beta, 0.0);

        // Sum over n
        int n = 1;
        bool converged = false;
        while (!converged) {
            // Compute coefficient
            double coef = pow(alpha * lambda_, n) / boost::math::tgamma(n + 1);

            std::vector<double> Tn;
            if (n == 1) {
                Tn = T1;
            } else {
                Tn = convolve(T1, Tn_list[n - 2]);
                Tn_list.push_back(Tn);
            }

            // Add to S_beta
            for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
                double term = coef * Tn[i_beta];
                S_beta[i_beta] += term;
            }

            // Check convergence
            double max_term = *max_element(Tn.begin(), Tn.end());
            if (coef * max_term < 1e-6) {
                converged = true;
            }

            n++;
            if (n > n_max_) {
                converged = true;
            }
        }

        // Multiply by exponential factor
        for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
            S_beta[i_beta] *= exp_factor;
        }

        // Store in S_alpha_beta_
        S_alpha_beta_[i_alpha] = S_beta;
    }
}
*/

/*
void SAlphaBeta::compute_S_alpha_beta()
{
    size_t n_alpha = alpha_grid_.size();
    size_t n_beta = beta_grid_.size();

    std::cout << "i_beta" << " "<< "n" << std::endl;

    // Create an interpolated rho(beta) function
    auto rho_beta_func = interpolate_rho_beta(beta_spectrum_grid_, rho_beta_);

    // Compute lambda
    auto P = [this, rho_beta_func](double beta) { return compute_P_beta(beta, rho_beta_func); };
    lambda_ = compute_lambda(P);
    compute_effective_temperature(P);

    auto T1_func = [this, rho_beta_func](double beta) { return compute_T1(beta, rho_beta_func); };

    // Initialize n_max matrix and alpha_max values
    std::vector<std::vector<int>> n_max_matrix(n_alpha, std::vector<int>(n_beta, 0));
    std::vector<double> alpha_max(n_beta, alpha_grid_.back());  // Default to last alpha

    // Compute T1(beta)
    std::vector<double> T1(n_beta, 0.0);
    for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
        T1[i_beta] = T1_func(beta_grid_[i_beta]);
    }

    std::vector<double> P1(n_beta, 0.0);
    for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
        P1[i_beta] = P(beta_grid_[i_beta]);
    }

    for (size_t i_alpha = 0; i_alpha < n_alpha; ++i_alpha) {
        double alpha = alpha_grid_[i_alpha];

        // Compute the first term: n=1 contribution
        std::vector<double> S_beta(n_beta, 0.0);
        double exp_factor = exp(-alpha * lambda_);  // e^(-alpha * lambda)
        for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
            // T1(beta) * e^(-alpha * lambda) * alpha * lambda
            S_beta[i_beta] = P1[i_beta] * alpha * std::exp( - beta_grid_[i_beta]/2 ); 
            if(i_beta == 0)  S_beta[i_beta] += 1;
        }

        // Sum over n, starting from n=2
        int n = 2;
        bool converged = false;
        while (!converged) {
            // Compute the nth term coefficient: (alpha^n / n!)
            double coef = pow(alpha, n) / boost::math::factorial<double>(n);

            // Compute Wn for the current n
            std::vector<double> Wn = (n == 2) ? convolve(P1, P1) : convolve(P1, Wn);

            // Check contributions and accumulate
            for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
                double term = coef * std::exp( - beta_grid_[i_beta]/2 ) * Wn[i_beta];  // Wn * e^(-beta/2)
                S_beta[i_beta] += term;

                // Check for convergence: if the term is less than 0.1% of the accumulated value
                if (term == 0 || std::abs(term) < 1e-2 * std::abs(S_beta[i_beta])) {
                    converged = true;
                    n_max_matrix[i_alpha][i_beta] = n;
                }
            }

            n++;
        }

        // Multiply by the overall exponential factor
        for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
            S_beta[i_beta] *= exp_factor;
        }

        // Store S(alpha, beta)
        S_alpha_beta_[i_alpha] = S_beta;
    }
}
*/

/*
std::vector<double> SAlphaBeta::convolve(const std::vector<double>& f, const std::vector<double>& g)
{
    // Perform convolution using numerical integration over the beta grid
    size_t N = beta_grid_.size();
    std::vector<double> result(N, 0.0);
    double beta_min_ = beta_grid_.front(); 
    double beta_max_ = beta_grid_.back(); 

    // Linear interpolation with symmetry handling
    auto interpolate_with_symmetry = [&](const std::vector<double>& vec, double beta_val) -> double {
        if (beta_val < beta_min_ || beta_val > beta_max_) {
            if (beta_val < 0 && -beta_val >= beta_min_ && -beta_val <= beta_max_) {
                // Use the symmetry T_n(beta) = e^{-beta} T_n(-beta)
                auto it = std::lower_bound(beta_grid_.begin(), beta_grid_.end(), -beta_val);
                size_t idx = std::distance(beta_grid_.begin(), it);

                if (idx == 0 || idx >= N) {
                    return vec[idx == 0 ? 0 : N - 1] * std::exp(beta_val);
                }

                // Linear interpolation for the positive counterpart -beta
                double beta_low = beta_grid_[idx - 1];
                double beta_high = beta_grid_[idx];
                double val_low = vec[idx - 1];
                double val_high = vec[idx];

                double weight = (-beta_val - beta_low) / (beta_high - beta_low);
                double interpolated_val = val_low + weight * (val_high - val_low);

                return interpolated_val * std::exp(beta_val);  // Apply the exponential factor e^{-\beta}
            }

            return 0.0; // Out of bounds
        }

        // Normal case: use linear interpolation as before
        auto it = std::lower_bound(beta_grid_.begin(), beta_grid_.end(), beta_val);
        size_t idx = std::distance(beta_grid_.begin(), it);

        if (idx == 0 || idx >= N) {
            return vec[idx == 0 ? 0 : N - 1]; // Boundary condition
        }

        double beta_low = beta_grid_[idx - 1];
        double beta_high = beta_grid_[idx];
        double val_low = vec[idx - 1];
        double val_high = vec[idx];

        double weight = (beta_val - beta_low) / (beta_high - beta_low);
        return val_low + weight * (val_high - val_low);
    };

    for (size_t i = 0; i < N; ++i) {
        double beta = beta_grid_[i];

        // Integrate over beta_prime
        double integral = boost::math::quadrature::trapezoidal(
            [&](double beta_prime) {
                double beta_minus_beta_prime = beta - beta_prime;

                // Interpolate f and g at beta_prime and beta - beta_prime
                double f_val = interpolate_with_symmetry(f, beta_prime);
                double g_val = interpolate_with_symmetry(g, beta_minus_beta_prime);

                return f_val * g_val;
            },
            beta_min_, beta_max_, 1e-6);

        result[i] = integral;
    }

    return result;
}
*/

void SAlphaBeta::compute_S_alpha_beta()
{
    auto rho_beta_func = interpolate_rho_beta(beta_spectrum_grid_, rho_beta_);
    auto P = [this, rho_beta_func](double beta) { return compute_P_beta(beta, rho_beta_func); };
    lambda_ = compute_lambda(P);
    compute_effective_temperature(P);

    for (size_t i_alpha = 0; i_alpha < alpha_grid_.size(); ++i_alpha) {
        double alpha = alpha_grid_[i_alpha];
        for (size_t i_beta = 0; i_beta < beta_grid_.size(); ++i_beta) {
            double beta = beta_grid_[i_beta];
            // Compute S(alpha, -beta)
            double alpha_term = 4 * w_s_ * alpha * T_s_bar_ / temperature_;
            double exp_term = std::exp( - (w_s_ * alpha - beta) * (w_s_ * alpha - beta) / alpha_term );

            double S_alpha_neg_beta = (1.0 / std::sqrt(4.0 * M_PI * alpha_term)) * exp_term;

            // Compute S(alpha, beta) using the relation S(alpha, beta) = e^(-beta) * S(alpha, -beta)
            S_alpha_beta_[i_alpha][i_beta] = std::exp( - beta ) * S_alpha_neg_beta;
        }
    }
}

/*
void SAlphaBeta::compute_S_alpha_beta()
{
    size_t n_alpha = alpha_grid_.size();
    size_t n_beta = beta_grid_.size();

    // Create an interpolated rho(beta) function
    auto rho_beta_func = interpolate_rho_beta(beta_spectrum_grid_, rho_beta_);

    // Compute lambda
    auto P = [this, rho_beta_func](double beta) { return compute_P_beta(beta, rho_beta_func); };
    lambda_ = compute_lambda(P);
    compute_effective_temperature(P);

    auto T1_func = [this, rho_beta_func](double beta) { return compute_T1(beta, rho_beta_func); };

    // Compute T1(beta) for all beta
    std::vector<double> T1(n_beta, 0.0);
    for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
        T1[i_beta] = T1_func(beta_grid_[i_beta]);
    }

    // For each alpha
    for (size_t i_alpha = 0; i_alpha < n_alpha; ++i_alpha) {
        double alpha = alpha_grid_[i_alpha];

        // Initialize S_beta for the current alpha
        std::vector<double> S_beta(n_beta, 0.0);

        // Exponential factor
        double exp_factor = std::exp(-alpha * lambda_);

        // Step 1: Add the n = 0 term (non-zero only at beta = 0)
        //S_beta[0] = 1;

        // Step 2: Add the n = 1 term (alpha * lambda * exp(-alpha * lambda) * T1(beta))
        for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
            S_beta[i_beta] += alpha * lambda_ * T1[i_beta];
        }

        // Step 3: Loop for n >= 2
        std::vector<double> Tn = convolve(T1, T1);  // T2 = convolve(T1, T1)

        for (int n = 2; n <= n_max_; ++n) {
            // Compute the coefficient (alpha * lambda)^n / n!
            double coef = std::pow(alpha * lambda_, n) / std::tgamma(n + 1);

            // Add the contribution of Tn to S_beta
            for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
                S_beta[i_beta] += coef * Tn[i_beta];
            }

            // Compute T_{n+1}(beta) by convolving T1(beta) with Tn(beta)
            Tn = convolve(T1, Tn);  // Tn gets updated to T_{n+1}(beta)
        }

        // Multiply by the overall exponential factor (already applied in step 1 and step 2)
        for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
            S_beta[i_beta] *= exp_factor;
        }

        // Store S_beta in S_alpha_beta_
        S_alpha_beta_[i_alpha] = S_beta;
    }
}
*/




// void SAlphaBeta::enforce_symmetry()
// {
//     size_t n_alpha = alpha_grid_.size();
//     size_t n_beta = beta_grid_.size();

//     for (size_t i_alpha = 0; i_alpha < n_alpha; ++i_alpha) {
//         for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
//             double beta = beta_grid_[i_beta];
//             double S_beta = S_alpha_beta_[i_alpha][i_beta];

//             // Find index of -beta
//             double neg_beta = -beta;
//             auto it_neg = lower_bound(beta_grid_.begin(), beta_grid_.end(), neg_beta);
//             if (it_neg == beta_grid_.end()) continue;
//             size_t i_beta_neg = distance(beta_grid_.begin(), it_neg);

//             double S_neg_beta = S_alpha_beta_[i_alpha][i_beta_neg];

//             // Enforce symmetry
//             double S_sym = 0.5 * (S_beta + exp(-beta) * S_neg_beta);
//             S_alpha_beta_[i_alpha][i_beta] = S_sym;
//             S_alpha_beta_[i_alpha][i_beta_neg] = exp(-beta) * S_sym;
//         }
//     }
// }

/*
void SAlphaBeta::compute_S_alpha_beta() 
{
    // Create an interpolated rho(beta) function
    auto rho_beta_func = interpolate_rho_beta(beta_spectrum_grid_, rho_beta_);
    
    // Compute lambda
    auto P = [this, rho_beta_func](double beta) { return compute_P_beta(beta, rho_beta_func); };

    std::ofstream output_file("P_values.txt");
    if (!output_file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Output P(beta) for a range of beta values
    for (double beta = 0.0; beta <= beta_spectrum_grid_.back(); beta += 0.01) {  // Adjust beta range as needed
        output_file << beta << "\t" << P(beta) << "\n";
    }

    output_file.close();

    lambda_ = compute_lambda(P);

    // Step 5: Compute effective temperature T_s_bar
    compute_effective_temperature(P);

    // Lambda function for T1
    auto T1_func = [this, rho_beta_func](double beta) { return compute_T1(beta, rho_beta_func); };

    output_debug_information();

    // Recursive function for Tn
    std::function<double(int, double)> Tn_func;

    Tn_func = [T1_func, &Tn_func](int n, double beta) {
        if (n == 1) {
            return T1_func(beta);
        } else {
            return compute_Tn(n, beta, T1_func, Tn_func);
        }
    };

    // Iterate over all alpha and beta in beta_grid_
    for (double alpha : beta_grid_) {
        for (double beta : beta_grid_) {
            double S_alpha_beta = compute_SAB(alpha, beta, n_max_, lambda_, Tn_func);
        }
    }
}
*/

void SAlphaBeta::determine_alpha_max(std::vector<std::vector<double>>& T_n)
{
    size_t N_beta = beta_grid_.size();
    size_t N_alpha = alpha_grid_.size();

    for (size_t beta_idx = 0; beta_idx < N_beta; ++beta_idx) {
        double S_sum = 0.0;
        bool converged = false;

        for (size_t alpha_idx = 0; alpha_idx < N_alpha; ++alpha_idx) {
            double alpha = alpha_grid_[alpha_idx];

            for (int n = 1; n <= n_max_; ++n) {
                if (n > 1) {
                    // Compute T_n(beta) using convolution
                    double convolution = 0.0;
                    for (size_t i = 0; i < N_beta; ++i) {
                        double beta_prime = beta_grid_[i];
                        double beta_diff = beta_grid_[beta_idx] - beta_prime;

                        // Find T1(beta - beta_prime)
                        double T1_value = 0.0;
                        if (beta_diff < beta_grid_.front() || beta_diff > beta_grid_.back()) {
                            T1_value = 0.0;
                        } else {
                            // Interpolate T1(beta_diff)
                            auto it = std::lower_bound(beta_grid_.begin(), beta_grid_.end(), beta_diff);
                            size_t idx = std::distance(beta_grid_.begin(), it);
                            if (idx == 0 || idx >= N_beta) {
                                T1_value = T1_beta_[idx];
                            } else {
                                double beta_low = beta_grid_[idx - 1];
                                double beta_high = beta_grid_[idx];
                                double T1_low = T1_beta_[idx - 1];
                                double T1_high = T1_beta_[idx];
                                double t = (beta_diff - beta_low) / (beta_high - beta_low);
                                T1_value = T1_low + t * (T1_high - T1_low);
                            }
                        }

                        double delta_beta = beta_grid_[1] - beta_grid_[0]; // Assuming uniform grid
                        convolution += T1_value * T_n[n - 1][i] * delta_beta;
                    }
                    T_n[n][beta_idx] = convolution;
                }

                // Compute Term
                double term = pow(alpha * lambda_, n) / tgamma(n + 1) * T_n[n][beta_idx];
                S_sum += term;

                // Check convergence
                if (fabs(term) < 0.001 * fabs(S_sum)) {
                    alpha_max_[beta_idx] = alpha;
                    converged = true;
                    break;
                }
            }
            if (converged) {
                break;
            }
        }
    }
}

/*
void SAlphaBeta::compute_S_alpha_beta_values(const std::vector<std::vector<double>>& T_n)
{
    size_t N_beta = beta_grid_.size();
    size_t N_alpha = alpha_grid_.size();

    for (size_t beta_idx = 0; beta_idx < N_beta; ++beta_idx) {
        double beta = beta_grid_[beta_idx];
        for (size_t alpha_idx = 0; alpha_idx < N_alpha; ++alpha_idx) {
            double alpha = alpha_grid_[alpha_idx];

            if (alpha <= alpha_max_[beta_idx]) {
                // Compute S(alpha, beta) using phonon expansion
                double S_sum = 0.0;
                for (int n = 1; n <= n_max_; ++n) {
                    double term = pow(alpha * lambda_, n) / tgamma(n + 1) * T_n[n][beta_idx];
                    S_sum += term;
                }
                double S_value = exp(-alpha * lambda_) * S_sum;
                S_alpha_beta_[beta_idx][alpha_idx] = S_value;
            } else {
                // Use SCT approximation
                double denom = w_s_ * alpha * T_s_bar_ / temperature_;
                if (denom <= 0.0) {
                    S_alpha_beta_[beta_idx][alpha_idx] = 0.0;
                    continue;
                }
                double sqrt_arg = 4.0 * M_PI * denom;
                double prefactor = 1.0 / sqrt(sqrt_arg);
                double exponent = -pow(w_s_ * alpha - beta, 2) / denom;
                double S_value = prefactor * exp(exponent);
                S_alpha_beta_[beta_idx][alpha_idx] = S_value;
            }
        }
    }
}
*/
 
//! @brief Returns the computed S(alpha, beta) matrix as a NumPy array.
//! @return The S(alpha, beta) as a pybind11 NumPy array.
//! 
// py::array_t<double> SAlphaBeta::get_S_alpha_beta()
// {
//     size_t n_alpha = alpha_grid_.size();
//     size_t n_beta = beta_grid_.size();

//     // Create NumPy array
//     py::array_t<double> S_array({n_alpha, n_beta});
//     auto S_buf = S_array.request();
//     double* S_ptr = static_cast<double*>(S_buf.ptr);

//     for (size_t i_alpha = 0; i_alpha < n_alpha; ++i_alpha) {
//         for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
//             S_ptr[i_alpha * n_beta + i_beta] = S_alpha_beta_[i_alpha][i_beta];
//         }
//     }

//     return S_array;
// }
py::array_t<double> SAlphaBeta::get_S_alpha_beta()
{
    // Get the dimensions of the S_alpha_beta_ matrix
    size_t n_alpha = S_alpha_beta_.size();
    size_t n_beta = S_alpha_beta_.empty() ? 0 : S_alpha_beta_[0].size();

    // Create a flat array to store the data in row-major order
    std::vector<double> flat_data;
    flat_data.reserve(n_alpha * n_beta);  // Reserve memory for the flattened array

    // Flatten the 2D S_alpha_beta_ into a 1D vector
    for (size_t i_alpha = 0; i_alpha < n_alpha; ++i_alpha) {
        for (size_t i_beta = 0; i_beta < n_beta; ++i_beta) {
            flat_data.push_back(S_alpha_beta_[i_alpha][i_beta]);
        }
    }

    // Create a py::array_t from the flat_data, with shape (n_alpha, n_beta)
    py::array_t<double> result({n_alpha, n_beta}, flat_data.data());

    return result;
}

void SAlphaBeta::output_debug_information()
{
    std::ofstream debug_output("debug_output.txt");
    if (!debug_output.is_open()) {
        std::cerr << "Failed to open debug_output.txt for writing.\n";
        return;
    }

    std::cout << "Lambda = " << lambda_ << std::endl;
    std::cout << "Eff. Temp. = " << T_s_bar_ << std::endl;

    // Print the header for the output file
    debug_output << "e\trho(e)\tbeta\trho(beta)\tP(beta)\tT1(beta)\n";

    // Lambda function for interpolated rho(beta)
    auto rho_beta_func = interpolate_rho_beta(beta_spectrum_grid_, rho_beta_);
    
    // Lambda function for P(beta)
    auto P_beta_func = [this, rho_beta_func](double beta) {
        return compute_P_beta(beta, rho_beta_func);
    };
    
    // Lambda function for T1(beta)
    auto T1_func = [this, P_beta_func](double beta) {
        return compute_T1(beta, P_beta_func);
    };

    size_t N_spectrum = energy_grid_.size();

    // Loop through all points in the spectrum
    for (size_t i = 0; i < N_spectrum; ++i) {
        double beta = beta_spectrum_grid_[i];

        // Calculate rho(beta), P(beta), and T1(beta) for the current beta
        double rho_beta = rho_beta_func(beta);
        double P_beta = P_beta_func(beta);
        double T1_beta = T1_func(beta);

        // Output all the computed values to the debug file
        debug_output << energy_grid_[i] << "\t" << rho_energy_[i] << "\t" << beta << "\t"
                     << rho_beta << "\t" << P_beta << "\t" << T1_beta << "\n";
    }

    debug_output.close();
    std::cout << "Debug information written to debug_output.txt.\n";
}