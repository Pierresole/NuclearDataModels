#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "ThermalScattering.h"

PYBIND11_MODULE(pyrat, m) {
    py::class_<SAlphaBeta>(m, "SAlphaBeta")
            .def(py::init<
                double,
                double, // Add w_s parameter
                const std::vector<double>&,
                const std::vector<double>&,
                int,
                const std::vector<double>&,
                const std::vector<double>&,
                const std::vector<std::tuple<double, double, std::string, double>>&>(),
                py::arg("temperature"),
                py::arg("w_s"), // Add w_s argument
                py::arg("alpha_grid"),
                py::arg("beta_grid"),
                py::arg("n_max"),
                py::arg("energy_grid"),
                py::arg("rho_energy"),
                py::arg("peaks") = std::vector<std::tuple<double, double, std::string, double>>())
            .def("compute_S_alpha_beta", &SAlphaBeta::compute_S_alpha_beta)
            .def("get_S_alpha_beta", &SAlphaBeta::get_S_alpha_beta)
            .def("get_alpha_grid", &SAlphaBeta::get_alpha_grid)
            .def("get_beta_grid", &SAlphaBeta::get_beta_grid)
            .def("get_lambda", &SAlphaBeta::get_lambda)
            .def("get_T_s_bar", &SAlphaBeta::get_T_s_bar)
            .def("get_alpha_max", &SAlphaBeta::get_alpha_max)
            .def("output_debug_information", &SAlphaBeta::output_debug_information);
}