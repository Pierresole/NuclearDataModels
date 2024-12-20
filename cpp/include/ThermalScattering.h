#ifndef SALPHABETA_H
#define SALPHABETA_H


#include <algorithm>
#include <cmath>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>
#include <stdexcept>
#include <string>
#include <tuple>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/special_functions/sinhc.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/constants/constants.hpp>

class SAlphaBeta
{
    double temperature_;                            // Temperature in Kelvin
    double w_s_;                                    // Continuous spectrum weight (provided)
    std::vector<double> beta_grid_;                 // Beta grid for calculations
    std::vector<double> beta_spectrum_grid_;        // Beta values derived from energy grid
    std::vector<double> alpha_grid_;                // Alpha grid
    std::vector<double> rho_beta_;                  // Phonon spectrum rho(beta)
    int n_max_;                                     // Max terms in phonon expansion

    double lambda_;                                 // Debye-Waller factor
    double T_s_bar_;                                // Effective temperature
    std::vector<double> alpha_max_;                 // Maximum alpha for each beta
    std::vector<std::vector<double>> S_alpha_beta_; // S(alpha, beta) values

    std::vector<double> P_beta_;
    std::vector<double> T1_beta_;

    // Input spectrum
    std::vector<double> energy_grid_;
    std::vector<double> rho_energy_;

    // Constants
    static constexpr double k_B = 8.617333262145e-5; // Boltzmann constant in eV/K

public:
    // Constructor
    SAlphaBeta(
        double temperature,
        double w_s, // New parameter
        const std::vector<double>& alpha_grid,
        const std::vector<double>& beta_grid,
        int n_max,
        const std::vector<double>& energy_grid,
        const std::vector<double>& rho_energy,
        const std::vector<std::tuple<double, double, std::string, double>>& peaks);

    // Getter methods
    const std::vector<double>& get_alpha_grid() const { return alpha_grid_; }
    const std::vector<double>& get_beta_grid() const { return beta_grid_; }
    double get_lambda() const { return lambda_; }
    double get_T_s_bar() const { return T_s_bar_; }
    const std::vector<double>& get_alpha_max() const { return alpha_max_; }
    py::array_t<double> get_S_alpha_beta();


    double compute_lambda(const std::function<double(double)>& P_beta_func);
    void compute_effective_temperature(const std::function<double(double)>& P_beta_func);
    double compute_P_beta(double beta, const std::function<double(double)>& rho_beta_func);
    double compute_T1(double beta, const std::function<double(double)>& rho_beta_func);
    void read_and_convert_spectrum();
    void determine_alpha_max(std::vector<std::vector<double>>& T_n);
    std::vector<double> convolve(const std::vector<double>&, const std::vector<double>&);
    void compute_S_alpha_beta();
    void output_debug_information();

    //Test in log space
    double interpolate_ln_T(const std::vector<double>& ln_T, double beta);
    std::vector<double> compute_ln_Tn(const std::vector<double>& ln_T1, const std::vector<double>& ln_Tn_minus1);

};

#endif // SALPHABETA_H
