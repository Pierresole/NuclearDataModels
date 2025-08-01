#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "CompoundSystem.h"
#include "SpinGroup.h"
#include "Channel.h"
#include "ParticlePair.h"
#include "FormalismParameters.h"
#include "Mlbw.h"  // Include the MLBW header
#include "ParameterSensitivity.h"

namespace py = pybind11;

// Helper class for simplified CompoundSystem construction
class CompoundSystemBuilder {
public:
    static CompoundSystem fromSpinGroups(const ParticlePair& entrancePP, const std::vector<SpinGroup>& spinGroups) {
        // Validate input
        if (spinGroups.empty()) {
            throw std::invalid_argument("No spin groups provided");
        }
        
        // Check that all spin groups are valid
        for (const auto& sg : spinGroups) {
            if (!sg.isValid()) {
                throw std::invalid_argument("Invalid spin group in construction");
            }
        }
        
        return CompoundSystem(entrancePP, spinGroups);
    }
    
    static CompoundSystem fromDict(const py::dict& data) {
        if (!data.contains("entrance_pair") || !data.contains("spin_groups")) {
            throw std::invalid_argument("Data dictionary must contain 'entrance_pair' and 'spin_groups'");
        }
        
        py::dict pp_data = data["entrance_pair"];
        ParticlePair entrancePP = createParticlePairFromDict(pp_data);
        
        std::vector<SpinGroup> spinGroups;
        for (const auto& sg_item : data["spin_groups"]) {
            py::dict sg_data = sg_item.cast<py::dict>();
            spinGroups.push_back(createSpinGroupFromDict(sg_data));
        }
        
        return fromSpinGroups(entrancePP, spinGroups);
    }
    
private:
    static ParticlePair createParticlePairFromDict(const py::dict& pp_data) {
        if (!pp_data.contains("mass1") || !pp_data.contains("mass2") || 
            !pp_data.contains("spin1") || !pp_data.contains("spin2") ||
            !pp_data.contains("Q") || !pp_data.contains("parity1") ||
            !pp_data.contains("parity2") || !pp_data.contains("MT")) {
            throw std::invalid_argument("Incomplete particle pair data");
        }
        
        return ParticlePair(
            pp_data["mass1"].cast<double>(),
            pp_data["mass2"].cast<double>(),
            pp_data["spin1"].cast<double>(),
            pp_data["spin2"].cast<double>(),
            pp_data["Q"].cast<double>(),
            pp_data["parity1"].cast<int>(),
            pp_data["parity2"].cast<int>(),
            pp_data["MT"].cast<int>()
        );
    }
    
    static SpinGroup createSpinGroupFromDict(const py::dict& sg_data) {
        if (!sg_data.contains("J") || !sg_data.contains("PJ")) {
            throw std::invalid_argument("Spin group must contain J and PJ");
        }
        
        SpinGroup sg(sg_data["J"].cast<double>(), sg_data["PJ"].cast<int>());
        
        if (sg_data.contains("channels")) {
            for (const auto& ch_item : sg_data["channels"]) {
                py::dict ch_data = ch_item.cast<py::dict>();
                sg.addChannel(createChannelFromDict(ch_data));
            }
        }
        
        if (sg_data.contains("resonances")) {
            for (const auto& res_item : sg_data["resonances"]) {
                py::dict res_data = res_item.cast<py::dict>();
                sg.addResonance(createResonanceFromDict(res_data));
            }
        }
        
        return sg;
    }
    
    static Channel createChannelFromDict(const py::dict& ch_data) {
        if (!ch_data.contains("particle_pair") || !ch_data.contains("l") || 
            !ch_data.contains("effective_radius") || !ch_data.contains("true_radius") ||
            !ch_data.contains("channel_spin")) {
            throw std::invalid_argument("Incomplete channel data");
        }
        
        ParticlePair pp = createParticlePairFromDict(ch_data["particle_pair"].cast<py::dict>());
        return Channel(
            pp,
            ch_data["l"].cast<unsigned int>(),
            ch_data["effective_radius"].cast<double>(),
            ch_data["true_radius"].cast<double>(),
            ch_data["channel_spin"].cast<double>()
        );
    }
    
    static Resonance createResonanceFromDict(const py::dict& res_data) {
        if (!res_data.contains("energy") || !res_data.contains("gammas")) {
            throw std::invalid_argument("Incomplete resonance data");
        }
        
        return Resonance(
            res_data["energy"].cast<double>(),
            res_data["gammas"].cast<std::vector<double>>()
        );
    }
};

// class CompoundWrapper {
// public:
//     CompoundSystem compound;

//     CompoundWrapper(const py::dict& data, double A) : compound(std::vector<SpinGroup>(), A) {
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

//             SpinGroup sg(AJ, PJ); // Use the second particle pair for now

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

//         compound = CrossSectionCalculator(spinGroups, A);
//     }

//     std::vector<double> computeCrossSections(const std::vector<double>& energies, const std::string& reactionChannel) {
//         return calculator.computeCrossSections(energies, reactionChannel);
//     }
// };

// Wrapper class to interface with Python
class FormalismParametersWrapper {
public:
    FormalismParameters data;

    FormalismParametersWrapper() : data() {}

    // Fill function to populate data from Python dict
    void fill(const py::dict& input_data) {
        data.SPI = input_data["SPI"].cast<double>();
        data.AP = input_data["AP"].cast<double>();
        data.LAD = input_data["LAD"].cast<int>();

        int LRF = input_data["LRF"].cast<int>();  // Read the flag to determine formalism

        if (LRF == 2) {  // MLBW Parameters
            MLBWParameters mlbwParams;
            mlbwParams.AWRI = input_data["AWRI"].cast<double>();
            mlbwParams.AP = input_data["AP"].cast<double>();
            auto llist = input_data["Llist"].cast<std::vector<py::dict>>();
            for (const auto& l_value : llist) {
                Llist llist_item;
                llist_item.APL = l_value["APL"].cast<double>();
                llist_item.L = l_value["L"].cast<int>();

                auto rplist = l_value["RPlist"].cast<py::dict>();
                RPlist rp_list;
                rp_list.ER = rplist["ER"].cast<std::vector<double>>();
                rp_list.AJ = rplist["AJ"].cast<std::vector<double>>();
                rp_list.GN = rplist["GN"].cast<std::vector<double>>();
                rp_list.GG = rplist["GG"].cast<std::vector<double>>();
                rp_list.GFA = rplist["GFA"].cast<std::vector<double>>();
                rp_list.GFB = rplist["GFB"].cast<std::vector<double>>();

                llist_item.resonanceParameters = rp_list;
                mlbwParams.lValues.push_back(llist_item);
            }
            data.parameters = mlbwParams;  // Set MLBW parameters
        } 
        else if (LRF == 7) {  // R-Matrix Parameters
            RMatrixParameters rmatrixParams;
            auto jpi_groups = input_data["JPiGroups"].cast<std::vector<py::dict>>();
            for (const auto& jpi_value : jpi_groups) {
                JPiList jpi_item;
                jpi_item.J = jpi_value["J"].cast<double>();
                jpi_item.PI = jpi_value["PI"].cast<double>();
                // Add other R-Matrix specific parameters if needed
                rmatrixParams.JPiGroups.push_back(jpi_item);
            }
            data.parameters = rmatrixParams;  // Set R-Matrix parameters
        }
    }
};

PYBIND11_MODULE(pyRMatrix, m) {

    py::class_<FormalismParametersWrapper>(m, "FormalismParametersWrapper")
        .def(py::init<>())
        .def("fill", &FormalismParametersWrapper::fill);  // Expose the fill method

    // Expose MLBW class
    py::class_<MLBW>(m, "MLBW")
        .def(py::init<const MLBWParameters&>())  // Constructor with parameters
        .def("radiative_capture_cross_section", &MLBW::radiative_capture_cross_section);  // Binding radiative capture cross-section

    // Optionally, expose structs if needed
    py::class_<RPlist>(m, "RPlist")
        .def(py::init<>())
        .def_readwrite("ER", &RPlist::ER)
        .def_readwrite("AJ", &RPlist::AJ)
        .def_readwrite("GN", &RPlist::GN)
        .def_readwrite("GG", &RPlist::GG)
        .def_readwrite("GFA", &RPlist::GFA)
        .def_readwrite("GFB", &RPlist::GFB);

    py::class_<Llist>(m, "Llist")
        .def(py::init<>())
        .def_readwrite("APL", &Llist::APL)
        .def_readwrite("L", &Llist::L)
        .def_readwrite("resonanceParameters", &Llist::resonanceParameters);

    py::class_<JPiList>(m, "JPiList")
        .def(py::init<>())
        .def_readwrite("J", &JPiList::J)
        .def_readwrite("PI", &JPiList::PI);

    py::class_<MLBWParameters>(m, "MLBWParameters")
        .def(py::init<>())
        .def_readwrite("AP", &MLBWParameters::AP)
        .def_readwrite("AWRI", &MLBWParameters::AWRI)
        .def_readwrite("lValues", &MLBWParameters::lValues);

    py::class_<RMatrixParameters>(m, "RMatrixParameters")
        .def(py::init<>())
        .def_readwrite("JPiGroups", &RMatrixParameters::JPiGroups);

    // Enhanced bindings for ParticlePair
    py::class_<ParticlePair>(m, "ParticlePair")
        .def(py::init<double, double, double, double, double, int, int, int>(),
             py::arg("mass1"), py::arg("mass2"), 
             py::arg("spin1"), py::arg("spin2"), 
             py::arg("Q"), 
             py::arg("parity1"), py::arg("parity2"), 
             py::arg("MT"))
        .def_static(
            "neutron_incident",
            &ParticlePair::neutronIncident,
            py::arg("mass2"),
            py::arg("spin2"),
            py::arg("Q"),
            py::arg("parity2"),
            py::arg("MT"),
            "Construct a ParticlePair with a neutron as the incident particle"
        )
        .def("mass1", &ParticlePair::mass1)
        .def("mass2", &ParticlePair::mass2)
        .def("spin1", &ParticlePair::spin1)
        .def("spin2", &ParticlePair::spin2)
        .def("MT", &ParticlePair::MT)
        .def("k2", &ParticlePair::k2)
        .def("reducedMass", &ParticlePair::reducedMass);

    // Enhanced bindings for Channel
    py::class_<Channel>(m, "Channel")
        .def(py::init<const ParticlePair&, unsigned int, double, double, double>(),
             py::arg("particle_pair"), py::arg("l"), 
             py::arg("effective_radius"), py::arg("true_radius"), 
             py::arg("channel_spin"))
        .def("computeChannelQuantities", &Channel::computeChannelQuantities)
        .def("isValid", &Channel::isValid)
        .def("L", &Channel::L)
        .def("getParticlePair", &Channel::getParticlePair, py::return_value_policy::reference_internal)
        // ...existing methods...
        ;
        
    // Enhanced bindings for Resonance
    py::class_<Resonance>(m, "Resonance")
        .def(py::init<double, const std::vector<double>&>(),
             py::arg("energy"), py::arg("gamma"))
        .def("getEnergy", &Resonance::getEnergy)
        .def("getGamma", &Resonance::getGamma);

    // Enhanced bindings for SpinGroup
    py::class_<SpinGroup>(m, "SpinGroup")
        .def(py::init<double, int>(), py::arg("J"), py::arg("PJ"))
        .def(py::init<double, int, const std::vector<Channel>&, const std::vector<Resonance>&>(),
             py::arg("J"), py::arg("PJ"), py::arg("channels"), py::arg("resonances"))
        .def("addChannel", &SpinGroup::addChannel)
        .def("addResonance", py::overload_cast<const Resonance&>(&SpinGroup::addResonance))
        .def("addResonance", py::overload_cast<double, const std::vector<double>&, bool>(&SpinGroup::addResonance),
             py::arg("energy"), py::arg("widths"), py::arg("isReduced") = true)
        .def("isValid", &SpinGroup::isValid)
        .def("channels", &SpinGroup::channels, py::return_value_policy::reference_internal)
        .def("getResonances", &SpinGroup::getResonances, py::return_value_policy::reference_internal)
        .def("getJ", &SpinGroup::getJ)
        .def("getPJ", &SpinGroup::getPJ)
        .def("crossSection", &SpinGroup::crossSection, py::arg("energy"), py::arg("entrance_particle_pair"))
        .def("elasticCrossSection", &SpinGroup::elasticCrossSection, py::arg("energy"), py::arg("entrance_particle_pair"))
        .def("captureCrossSection", &SpinGroup::captureCrossSection, py::arg("energy"), py::arg("entrance_particle_pair"))
        .def("fissionCrossSection", &SpinGroup::fissionCrossSection, py::arg("energy"), py::arg("entrance_particle_pair"))
        .def("totalCrossSection", &SpinGroup::totalCrossSection, py::arg("energy"), py::arg("entrance_particle_pair"))
        // Overloaded methods for energy lists
        .def("crossSection", [](const SpinGroup& self, const std::vector<double>& energies, const ParticlePair& pp) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.crossSection(E, pp));
            }
            return results;
        }, py::arg("energies"), py::arg("entrance_particle_pair"))
        .def("elasticCrossSection", [](const SpinGroup& self, const std::vector<double>& energies, const ParticlePair& pp) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.elasticCrossSection(E, pp));
            }
            return results;
        }, py::arg("energies"), py::arg("entrance_particle_pair"))
        .def("captureCrossSection", [](const SpinGroup& self, const std::vector<double>& energies, const ParticlePair& pp) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.captureCrossSection(E, pp));
            }
            return results;
        }, py::arg("energies"), py::arg("entrance_particle_pair"))
        .def("fissionCrossSection", [](const SpinGroup& self, const std::vector<double>& energies, const ParticlePair& pp) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.fissionCrossSection(E, pp));
            }
            return results;
        }, py::arg("energies"), py::arg("entrance_particle_pair"))
        .def("totalCrossSection", [](const SpinGroup& self, const std::vector<double>& energies, const ParticlePair& pp) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.totalCrossSection(E, pp));
            }
            return results;
        }, py::arg("energies"), py::arg("entrance_particle_pair"));

    // Enhanced bindings for CompoundSystem with builder
    py::class_<CompoundSystem>(m, "CompoundSystem")
        .def(py::init<const ParticlePair&, const std::vector<SpinGroup>&>(), 
             py::arg("entrance_particle_pair"), py::arg("spin_groups"))
        .def("addSpinGroup", &CompoundSystem::addSpinGroup)
        .def("printSpinGroupInfo", &CompoundSystem::printSpinGroupInfo)
        .def("spinGroups", [](const CompoundSystem& self) {
            auto spinGroupsView = self.spinGroups();
            return std::vector<SpinGroup>(spinGroupsView.begin(), spinGroupsView.end());
        }, "Retrieve spin groups as a Python-compatible list")
        .def("crossSection", &CompoundSystem::crossSection, py::arg("energy"))
        .def("elasticCrossSection", &CompoundSystem::elasticCrossSection, py::arg("energy"))
        .def("captureCrossSection", &CompoundSystem::captureCrossSection, py::arg("energy"))
        .def("fissionCrossSection", &CompoundSystem::fissionCrossSection, py::arg("energy"))
        .def("totalCrossSection", &CompoundSystem::totalCrossSection, py::arg("energy"))
        .def("spinGroupElasticCrossSection", &CompoundSystem::spinGroupElasticCrossSection, 
             py::arg("spin_group_index"), py::arg("energy"))
        .def("spinGroupCaptureCrossSection", &CompoundSystem::spinGroupCaptureCrossSection, 
             py::arg("spin_group_index"), py::arg("energy"))
        .def("spinGroupFissionCrossSection", &CompoundSystem::spinGroupFissionCrossSection, 
             py::arg("spin_group_index"), py::arg("energy"))
        .def("spinGroupTotalCrossSection", &CompoundSystem::spinGroupTotalCrossSection, 
             py::arg("spin_group_index"), py::arg("energy"))
        // Overloaded methods for energy lists
        .def("crossSection", [](const CompoundSystem& self, const std::vector<double>& energies) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.crossSection(E));
            }
            return results;
        }, py::arg("energies"))
        .def("elasticCrossSection", [](const CompoundSystem& self, const std::vector<double>& energies) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.elasticCrossSection(E));
            }
            return results;
        }, py::arg("energies"))
        .def("captureCrossSection", [](const CompoundSystem& self, const std::vector<double>& energies) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.captureCrossSection(E));
            }
            return results;
        }, py::arg("energies"))
        .def("fissionCrossSection", [](const CompoundSystem& self, const std::vector<double>& energies) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.fissionCrossSection(E));
            }
            return results;
        }, py::arg("energies"))
        .def("totalCrossSection", [](const CompoundSystem& self, const std::vector<double>& energies) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.totalCrossSection(E));
            }
            return results;
        }, py::arg("energies"))
        .def("spinGroupElasticCrossSection", [](const CompoundSystem& self, size_t spinGroupIndex, const std::vector<double>& energies) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.spinGroupElasticCrossSection(spinGroupIndex, E));
            }
            return results;
        }, py::arg("spin_group_index"), py::arg("energies"))
        .def("spinGroupCaptureCrossSection", [](const CompoundSystem& self, size_t spinGroupIndex, const std::vector<double>& energies) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.spinGroupCaptureCrossSection(spinGroupIndex, E));
            }
            return results;
        }, py::arg("spin_group_index"), py::arg("energies"))
        .def("spinGroupFissionCrossSection", [](const CompoundSystem& self, size_t spinGroupIndex, const std::vector<double>& energies) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.spinGroupFissionCrossSection(spinGroupIndex, E));
            }
            return results;
        }, py::arg("spin_group_index"), py::arg("energies"))
        .def("spinGroupTotalCrossSection", [](const CompoundSystem& self, size_t spinGroupIndex, const std::vector<double>& energies) {
            std::vector<double> results;
            results.reserve(energies.size());
            for (double E : energies) {
                results.push_back(self.spinGroupTotalCrossSection(spinGroupIndex, E));
            }
            return results;
        }, py::arg("spin_group_index"), py::arg("energies"))
        .def("entranceParticlePair", &CompoundSystem::entranceParticlePair)
        .def("getSpinGroup", &CompoundSystem::getSpinGroup, py::return_value_policy::reference_internal);

    // py::class_<CompoundWrapper>(m, "CrossSectionCalculator")
    //     .def(py::init<const py::dict&, double>(), py::arg("data"), py::arg("A"))
    //     .def("computeCrossSections", &CompoundWrapper::crossSections,
    //          py::arg("energies"), py::arg("reactionChannel"));

        
    // Add CompoundSystemBuilder for easy construction
    py::class_<CompoundSystemBuilder>(m, "CompoundSystemBuilder")
        .def_static("from_spin_groups", &CompoundSystemBuilder::fromSpinGroups, 
                    py::arg("entrance_particle_pair"), py::arg("spin_groups"))
        .def_static("from_dict", &CompoundSystemBuilder::fromDict, 
                    py::arg("data"));
    
    // Bindings for ParameterInfo
    py::class_<ParameterInfo>(m, "ParameterInfo")
        .def_readwrite("type", &ParameterInfo::type)
        .def_readwrite("spinGroupIndex", &ParameterInfo::spinGroupIndex)
        .def_readwrite("resonanceIndex", &ParameterInfo::resonanceIndex)
        .def_readwrite("channelIndex", &ParameterInfo::channelIndex)
        .def_readwrite("particlePairIndex", &ParameterInfo::particlePairIndex)
        .def_readwrite("description", &ParameterInfo::description)
        .def_readwrite("nominalValue", &ParameterInfo::nominalValue);
    
    // Enum for ParameterInfo::Type
    py::enum_<ParameterInfo::Type>(m, "ParameterType")
        .value("RESONANCE_ENERGY", ParameterInfo::Type::RESONANCE_ENERGY)
        .value("RESONANCE_GAMMA", ParameterInfo::Type::RESONANCE_GAMMA)
        .value("CHANNEL_RADIUS", ParameterInfo::Type::CHANNEL_RADIUS)
        .value("PARTICLE_MASS", ParameterInfo::Type::PARTICLE_MASS)
        .value("PARTICLE_SPIN", ParameterInfo::Type::PARTICLE_SPIN);
    
    // Bindings for ParameterSensitivity
    py::class_<ParameterSensitivity>(m, "ParameterSensitivity")
        .def(py::init<const CompoundSystem&, double>(), 
             py::arg("system"), py::arg("perturbation") = 1e-6)
        .def("extractParameters", &ParameterSensitivity::extractParameters)
        .def("computeDerivative", &ParameterSensitivity::computeDerivative,
             py::arg("parameter_index"), py::arg("energy"), py::arg("cross_section_func"), py::arg("perturbation") = -1.0)
        .def("totalCrossSectionDerivatives", &ParameterSensitivity::totalCrossSectionDerivatives,
             py::arg("energy"), py::arg("perturbation") = -1.0)
        .def("elasticCrossSectionDerivatives", &ParameterSensitivity::elasticCrossSectionDerivatives,
             py::arg("energy"), py::arg("perturbation") = -1.0)
        .def("captureCrossSectionDerivatives", &ParameterSensitivity::captureCrossSectionDerivatives,
             py::arg("energy"), py::arg("perturbation") = -1.0)
        .def("fissionCrossSectionDerivatives", &ParameterSensitivity::fissionCrossSectionDerivatives,
             py::arg("energy"), py::arg("perturbation") = -1.0)
        .def("getParameters", &ParameterSensitivity::getParameters, py::return_value_policy::reference_internal)
        .def("getParameterCount", &ParameterSensitivity::getParameterCount)
        .def("printParameterInfo", &ParameterSensitivity::printParameterInfo)
        .def("computeSensitivityMatrix", &ParameterSensitivity::computeSensitivityMatrix,
             py::arg("energies"), py::arg("cross_section_types") = std::vector<std::string>{"total", "elastic", "capture", "fission"});
}