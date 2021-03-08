/*
 * MIT License
 *
 * Copyright (c) 2018 Guillaume St-Onge
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <PowerlawGroupSIS.hpp>
#include <MarginalInfectionProbability.hpp>
#include <Prevalence.hpp>

using namespace std;
using namespace schon;

namespace py = pybind11;


PYBIND11_MODULE(_schon, m)
{
    /* ========================
     * Base class GroupSIS
     * ========================*/

    py::class_<GroupSIS>(m, "GroupSIS")

        .def(py::init<EdgeList, double,
                const function<double(size_t,size_t)>&,
                const pair<double,double>&>(), R"pbdoc(
            Default constructor of the class GroupSIS.

            Args:
               edge_list: Edge list for the network structure.
               recovery_rate: Double for the recovery rate
               infection_rate: Function for the recovery rate
               rate_bounds: Rate lower and upper bounds.
            )pbdoc", py::arg("edge_list"),
                py::arg("recovery_rate"),
                py::arg("infection_rate"),
                py::arg("rate_bounds"))

        .def("size", &GroupSIS::size, R"pbdoc(
            Returns the number of nodes.
            )pbdoc")

        .def("get_state_vector", &GroupSIS::get_node_state_vector, R"pbdoc(
            Returns the vector of state for each node.
            )pbdoc")

        .def("get_current_time", &GroupSIS::get_current_time, R"pbdoc(
            Returns the current time for the process.
            )pbdoc")

        .def("get_lifetime", &GroupSIS::get_lifetime, R"pbdoc(
            Returns the lifetime for the current state.
            )pbdoc")

        .def("get_number_of_infected_nodes",
                &GroupSIS::get_number_of_infected_nodes, R"pbdoc(
            Returns the number of infected nodes.
            )pbdoc")

        .def("infect_fraction", &GroupSIS::infect_fraction, R"pbdoc(
            Infect a fraction of the nodes.

            Args:
               fraction: Fraction to be infected.
            )pbdoc", py::arg("fraction"))

        .def("infect_node_set", &GroupSIS::infect_fraction, R"pbdoc(
            Infect the nodes in the node set.

            Args:
               node_set: Set of nodes to infect.
            )pbdoc", py::arg("node_set"))


        .def("evolve", &GroupSIS::evolve,
                R"pbdoc(
            Let the system evolve over a period of time.

            Args:
               period: Time period of the evolution.
               decorrelation_time (optional): Time period for decorrelation.
               measure: Bool, if true measure the system after decorrelation.
               quasistationary: Bool, if true quasistationary state.
            )pbdoc", py::arg("period"), py::arg("decorrelation_time")=1,
                py::arg("measure")=false, py::arg("quasistationary")=false)

        .def("clear", &GroupSIS::clear, R"pbdoc(
            Recover all nodes.
            )pbdoc")

        .def("reset", &GroupSIS::reset, R"pbdoc(
            Reset time and system, with all susceptible nodes, remove history.
            )pbdoc")

        .def("initialize_history", &GroupSIS::initialize_history,
                R"pbdoc(
            Fill the history with the current state.

            Args:
               number_of_states: Number of state copies.
            )pbdoc", py::arg("number_of_states")=100)

        .def("seed", &GroupSIS::seed,
                R"pbdoc(
            Seed the RNG.

            Args:
               seed: seed for the RNG.
            )pbdoc", py::arg("seed"))

        .def("get_measure_vector",
                &MeasurableContagionProcess::get_measure_vector, R"pbdoc(
            Returns a vector of pointer to different measures.
            )pbdoc")

        .def("measure_marginal_infection_probability",
        &MeasurableContagionProcess::measure_marginal_infection_probability,
            R"pbdoc(
            Add the measure of marginal probability of infection to the vector
            of measures to be performed during simulations.
            )pbdoc")

        .def("measure_prevalence",
        &MeasurableContagionProcess::measure_prevalence,
            R"pbdoc(
            Add the measure of prevalence to the vector
            of measures to be performed during simulations.
            )pbdoc");

    /* ===============================
     * Class deriving from GroupSIS
     * ===============================*/

    py::class_<PowerlawGroupSIS, GroupSIS>(m, "PowerlawGroupSIS")

        .def(py::init<EdgeList, double, double, double,
                const pair<double,double>&>(), R"pbdoc(
            Default constructor of the class PowerlawGroupSIS.

            Args:
               edge_list: Edge list for the network structure.
               scale_recovery: Recovery rate for infected nodes
               scale_infection: Infection rate factor.
               shape_infection: Power-law exponent for infection.
               rate_bounds: Rate lower and upper bounds.
            )pbdoc", py::arg("edge_list"),
                py::arg("scale_recovery"),
                py::arg("scale_infection"),
                py::arg("shape_infection"),
                py::arg("rate_bounds"));


    /* ================
     * Measure classes
     * ================*/

    py::class_<MarginalInfectionProbability,
        shared_ptr<MarginalInfectionProbability>>(m,
                "MarginalInfectionProbability")

        .def(py::init<std::size_t>(), R"pbdoc(
            Default constructor of the class MarginalInfectionProbability.

            Args:
               network_size: Number of nodes in the network.
            )pbdoc", py::arg("network_size"))

        .def("get_name", &MarginalInfectionProbability::get_name, R"pbdoc(
            Returns the name of the measure.
            )pbdoc")

        .def("get_result", &MarginalInfectionProbability::get_result, R"pbdoc(
            Returns the result associated to the measure.
            )pbdoc");

    py::class_<Prevalence,shared_ptr<Prevalence>>(m,
                "Prevalence")

        .def(py::init<std::size_t>(), R"pbdoc(
            Default constructor of the class Prevalence.

            Args:
               network_size: Number of nodes in the network.
            )pbdoc", py::arg("network_size"))

        .def("get_name", &Prevalence::get_name, R"pbdoc(
            Returns the name of the measure.
            )pbdoc")

        .def("get_result", &Prevalence::get_result, R"pbdoc(
            Returns the result associated to the measure.
            )pbdoc");
}
