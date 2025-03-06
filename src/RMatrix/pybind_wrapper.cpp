// // pybind_wrapper.cpp
// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
// #include "FormalismParameters.h"
// #include "ThermalScattering.h"
// #include "Mlbw.h"  // Include the MLBW header
// #include "CrossSectionCalculator.h"
// #include "ParticlePair.h"
// #include "Channel.h"
// #include "SpinGroup.h"
// #include "CompoundSystem.h"

// namespace py = pybind11;

// // Wrapper class to interface with Python
// class CrossSectionCalculatorWrapper {
// public:
//     CrossSectionCalculator calculator;

//     CrossSectionCalculatorWrapper(const py::dict& data, double A) : calculator(std::vector<SpinGroup>(), A) {
//         // Parse the data dict and construct the spinGroups vector
//         std::vector<SpinGroup> spinGroups;

//         // Get particle pairs
//         std::vector<ParticlePair> particlePairs;
//         for (const auto& pp_data : data["particle_pairs"].cast<py::list>()) {
//             py::dict pp_dict = pp_data.cast<py::dict>();
//             double IA = pp_dict["IA"].cast<double>();
//             double IB = pp_dict["IB"].cast<double>();
//             double MA = pp_dict["MA"].cast<double>();
//             double MB = pp_dict["MB"].cast<double>();
//             int MT = pp_dict["MT"].cast<int>();
//             int PA = pp_dict["PA"].cast<int>();
//             int PB = pp_dict["PB"].cast<int>();
//             int PNT = pp_dict["PNT"].cast<int>();
//             double Q = pp_dict["Q"].cast<double>();
//             int SHF = pp_dict["SHF"].cast<int>();
//             int ZA = pp_dict["ZA"].cast<int>();
//             int ZB = pp_dict["ZB"].cast<int>();

//             std::string reactionID;
//             // Determine reaction ID based on MT
//             if (MT == 2) {
//                 reactionID = "(n,n)";
//             } else if (MT == 102) {
//                 reactionID = "(n,g)";
//             } else {
//                 reactionID = "(n,other)";
//             }

//             ParticlePair pp(1.0, MB, IA, IB, Q, PA, PB, MT); // massParticleA_ = 1.0 (neutron)
//             particlePairs.push_back(pp);
//         }

//         // Get spin groups
//         for (const auto& sg_data : data["spin_groups"].cast<py::list>()) {
//             py::dict sg_dict = sg_data.cast<py::dict>();
//             double AJ = sg_dict["AJ"].cast<double>();
//             int PJ = sg_dict["PJ"].cast<int>();

//             SpinGroup sg(AJ, PJ, particlePairs[1]); // Use the second particle pair for now

//             // Get channels
//             for (const auto& ch_data : sg_dict["channels"].cast<py::list>()) {
//                 py::dict ch_dict = ch_data.cast<py::dict>();
//                 int PPI = ch_dict["PPI"].cast<int>();
//                 double SCH = ch_dict["SCH"].cast<double>();
//                 double APE = ch_dict["APE"].cast<double>();
//                 unsigned int L = ch_dict["L"].cast<unsigned int>();

//                 // Get the corresponding ParticlePair
//                 const ParticlePair& pp = particlePairs[PPI];

//                 // ChannelRadius radius(APE); // Let it be a double for now
//                 // Channel channel(pp, L, APE);

//                 // sg.addChannel(channel);
//             }

//             // Get resonances
//             for (const auto& res_data : sg_dict["resonances"].cast<py::list>()) {
//                 py::dict res_dict = res_data.cast<py::dict>();
//                 double ER = res_dict["ER"].cast<double>();
//                 std::vector<double> GAM = res_dict["GAM"].cast<std::vector<double>>();

//                 Resonance resonance(ER, GAM);
//                 sg.addResonance(resonance);
//             }

//             spinGroups.push_back(sg);
//         }

//         calculator = CrossSectionCalculator(spinGroups, A);
//     }

//     std::vector<double> computeCrossSections(const std::vector<double>& energies, const std::string& reactionChannel) {
//         return calculator.computeCrossSections(energies, reactionChannel);
//     }
// };

// // Wrapper class to interface with Python
// class FormalismParametersWrapper {
// public:
//     FormalismParameters data;

//     FormalismParametersWrapper() : data() {}

//     // Fill function to populate data from Python dict
//     void fill(const py::dict& input_data) {
//         data.SPI = input_data["SPI"].cast<double>();
//         data.AP = input_data["AP"].cast<double>();
//         data.LAD = input_data["LAD"].cast<int>();

//         int LRF = input_data["LRF"].cast<int>();  // Read the flag to determine formalism

//         if (LRF == 2) {  // MLBW Parameters
//             MLBWParameters mlbwParams;
//             mlbwParams.AWRI = input_data["AWRI"].cast<double>();
//             mlbwParams.AP = input_data["AP"].cast<double>();
//             auto llist = input_data["Llist"].cast<std::vector<py::dict>>();
//             for (const auto& l_value : llist) {
//                 Llist llist_item;
//                 llist_item.APL = l_value["APL"].cast<double>();
//                 llist_item.L = l_value["L"].cast<int>();

//                 auto rplist = l_value["RPlist"].cast<py::dict>();
//                 RPlist rp_list;
//                 rp_list.ER = rplist["ER"].cast<std::vector<double>>();
//                 rp_list.AJ = rplist["AJ"].cast<std::vector<double>>();
//                 rp_list.GN = rplist["GN"].cast<std::vector<double>>();
//                 rp_list.GG = rplist["GG"].cast<std::vector<double>>();
//                 rp_list.GFA = rplist["GFA"].cast<std::vector<double>>();
//                 rp_list.GFB = rplist["GFB"].cast<std::vector<double>>();

//                 llist_item.resonanceParameters = rp_list;
//                 mlbwParams.lValues.push_back(llist_item);
//             }
//             data.parameters = mlbwParams;  // Set MLBW parameters
//         } 
//         else if (LRF == 7) {  // R-Matrix Parameters
//             RMatrixParameters rmatrixParams;
//             auto jpi_groups = input_data["JPiGroups"].cast<std::vector<py::dict>>();
//             for (const auto& jpi_value : jpi_groups) {
//                 JPiList jpi_item;
//                 jpi_item.J = jpi_value["J"].cast<double>();
//                 jpi_item.PI = jpi_value["PI"].cast<double>();
//                 // Add other R-Matrix specific parameters if needed
//                 rmatrixParams.JPiGroups.push_back(jpi_item);
//             }
//             data.parameters = rmatrixParams;  // Set R-Matrix parameters
//         }
//     }
// };

// PYBIND11_MODULE(pyrat, m) {
//     py::class_<FormalismParametersWrapper>(m, "FormalismParametersWrapper")
//         .def(py::init<>())
//         .def("fill", &FormalismParametersWrapper::fill);  // Expose the fill method

//     // Expose MLBW class
//     py::class_<MLBW>(m, "MLBW")
//         .def(py::init<const MLBWParameters&>())  // Constructor with parameters
//         .def("radiative_capture_cross_section", &MLBW::radiative_capture_cross_section);  // Binding radiative capture cross-section

//     // Optionally, expose structs if needed
//     py::class_<RPlist>(m, "RPlist")
//         .def(py::init<>())
//         .def_readwrite("ER", &RPlist::ER)
//         .def_readwrite("AJ", &RPlist::AJ)
//         .def_readwrite("GN", &RPlist::GN)
//         .def_readwrite("GG", &RPlist::GG)
//         .def_readwrite("GFA", &RPlist::GFA)
//         .def_readwrite("GFB", &RPlist::GFB);

//     py::class_<Llist>(m, "Llist")
//         .def(py::init<>())
//         .def_readwrite("APL", &Llist::APL)
//         .def_readwrite("L", &Llist::L)
//         .def_readwrite("resonanceParameters", &Llist::resonanceParameters);

//     py::class_<JPiList>(m, "JPiList")
//         .def(py::init<>())
//         .def_readwrite("J", &JPiList::J)
//         .def_readwrite("PI", &JPiList::PI);

//     py::class_<MLBWParameters>(m, "MLBWParameters")
//         .def(py::init<>())
//         .def_readwrite("AP", &MLBWParameters::AP)
//         .def_readwrite("AWRI", &MLBWParameters::AWRI)
//         .def_readwrite("lValues", &MLBWParameters::lValues);

//     py::class_<RMatrixParameters>(m, "RMatrixParameters")
//         .def(py::init<>())
//         .def_readwrite("JPiGroups", &RMatrixParameters::JPiGroups);

//     py::class_<SAlphaBeta>(m, "SAlphaBeta")
//         .def(py::init<
//             double,
//             double, // Add w_s parameter
//             const std::vector<double>&,
//             const std::vector<double>&,
//             int,
//             const std::vector<double>&,
//             const std::vector<double>&,
//             const std::vector<std::tuple<double, double, std::string, double>>&>(),
//             py::arg("temperature"),
//             py::arg("w_s"), // Add w_s argument
//             py::arg("alpha_grid"),
//             py::arg("beta_grid"),
//             py::arg("n_max"),
//             py::arg("energy_grid"),
//             py::arg("rho_energy"),
//             py::arg("peaks") = std::vector<std::tuple<double, double, std::string, double>>())
//         .def("compute_S_alpha_beta", &SAlphaBeta::compute_S_alpha_beta)
//         .def("get_S_alpha_beta", &SAlphaBeta::get_S_alpha_beta)
//         .def("get_alpha_grid", &SAlphaBeta::get_alpha_grid)
//         .def("get_beta_grid", &SAlphaBeta::get_beta_grid)
//         .def("get_lambda", &SAlphaBeta::get_lambda)
//         .def("get_T_s_bar", &SAlphaBeta::get_T_s_bar)
//         .def("get_alpha_max", &SAlphaBeta::get_alpha_max)
//         .def("output_debug_information", &SAlphaBeta::output_debug_information);

//     py::class_<CrossSectionCalculatorWrapper>(m, "CrossSectionCalculator")
//         .def(py::init<const py::dict&, double>(), py::arg("data"), py::arg("A"))
//         .def("computeCrossSections", &CrossSectionCalculatorWrapper::computeCrossSections,
//              py::arg("energies"), py::arg("reactionChannel"));

//     // Expose other classes if needed
//     py::class_<ParticlePair>(m, "ParticlePair")
//         .def(py::init<double, double, double, double, double, int, int, int>(),
//              py::arg("massA"), py::arg("massB"), py::arg("spinA"), py::arg("spinB"), py::arg("QI"),
//              py::arg("parityA"), py::arg("parityB"), py::arg("MT"))
//         .def("getReactionID", &ParticlePair::getReactionID);

//     py::class_<Channel>(m, "Channel")
//         .def(py::init<const ParticlePair&, unsigned int, double, double, double>(),
//              py::arg("particlePair"), py::arg("l"), py::arg("effectiveChannelRadius"), py::arg("trueChannelRadius"), py::arg("channelSpin"))
//         .def("computeChannelQuantities", &Channel::computeChannelQuantities, py::arg("E"))
//         .def("PhaseShift", &Channel::PhaseShift, py::arg("energy"), py::arg("l"))
//         .def("PenetrationFactor", &Channel::PenetrationFactor, py::arg("energy"), py::arg("l"))
//         .def("ShiftFactor", &Channel::ShiftFactor, py::arg("energy"), py::arg("l"))
//         .def("getWaveNumber", &Channel::getWaveNumber)
//         .def("getWaveNumberSquared", &Channel::getWaveNumberSquared)
//         .def("getL", &Channel::getL)
//         .def("getParticlePair", &Channel::getParticlePair, py::return_value_policy::reference_internal);

//     py::class_<Resonance>(m, "Resonance")
//         .def(py::init<double, const std::vector<double>&>(),
//              py::arg("energy"), py::arg("gamma"))
//         .def("getEnergy", &Resonance::getEnergy)
//         .def("getGamma", &Resonance::getGamma);

//     py::class_<SpinGroup>(m, "SpinGroup")
//         .def(py::init<double, int>(), py::arg("J"), py::arg("PJ"))
//         .def(py::init<double, int, const std::vector<Channel>&, const std::vector<Resonance>&>(),
//              py::arg("J"), py::arg("PJ"), py::arg("channels"), py::arg("resonances"))
//         .def("addChannel", &SpinGroup::addChannel)
//         .def("addResonance", &SpinGroup::addResonance)
//         .def("channels", &SpinGroup::channels, py::return_value_policy::reference_internal)
//         .def("getResonances", &SpinGroup::getResonances, py::return_value_policy::reference_internal)
//         .def("getJ", &SpinGroup::getJ)
//         .def("getPJ", &SpinGroup::getPJ);
    
//     py::class_<CompoundSystem>(m, "CompoundSystem")
//         .def(py::init<const ParticlePair&>())
//         .def("addSpinGroup", &CompoundSystem::addSpinGroup)
//         .def("spinGroups", &CompoundSystem::spinGroups, py::return_value_policy::reference_internal)
//         .def("computeTotalCrossSection", &CompoundSystem::computeTotalCrossSection);
// }
