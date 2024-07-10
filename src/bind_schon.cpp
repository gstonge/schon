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
#include <BaseContagion.hpp>
#include <ContinuousSIS.hpp>
#include <ContinuousSIR.hpp>
#include <DiscreteSIS.hpp>
#include <HeterogeneousExposure.hpp>
#include <MarginalInfectionProbability.hpp>
#include <Prevalence.hpp>
#include <InfectiousSet.hpp>
#include <Time.hpp>

using namespace std;
using namespace schon;

namespace py = pybind11;


PYBIND11_MODULE(_schon, m)
{
    /* ===========
     * Base class
     * ===========*/

    py::class_<BaseContagion>(m, "BaseContagion")

        .def(py::init<EdgeList>(), R"pbdoc(
            Default constructor of the class BaseContagion.

            Args:
               edge_list: Edge list for the network structure.
            )pbdoc", py::arg("edge_list"))

        .def("size", &BaseContagion::size, R"pbdoc(
            Returns the number of nodes.
            )pbdoc")

        .def("get_state_vector", &BaseContagion::get_node_state_vector, R"pbdoc(
            Returns the vector of state for each node.
            )pbdoc")

        .def("get_current_time", &BaseContagion::get_current_time, R"pbdoc(
            Returns the current time for the process.
            )pbdoc")

        .def("get_number_of_infected_nodes",
                &BaseContagion::get_number_of_infected_nodes, R"pbdoc(
            Returns the number of infected nodes.
            )pbdoc")

        .def("infect_fraction", &BaseContagion::infect_fraction, R"pbdoc(
            Infect a fraction of the nodes.

            Args:
               fraction: Fraction to be infected.
            )pbdoc", py::arg("fraction"))

        .def("infect_node_set", &BaseContagion::infect_node_set, R"pbdoc(
            Infect the nodes in the node set.

            Args:
               node_set: Set of nodes to infect.
            )pbdoc", py::arg("node_set"))

        .def("clear", &BaseContagion::clear, R"pbdoc(
            Recover all nodes.
            )pbdoc")

        .def("reset", &BaseContagion::reset, R"pbdoc(
            Reset time and system, with all susceptible nodes, remove history.
            )pbdoc")

        .def("initialize_history", &BaseContagion::initialize_history,
                R"pbdoc(
            Fill the history with the current state.

            Args:
               number_of_states: Number of state copies.
            )pbdoc", py::arg("number_of_states")=100)

        .def("seed", &BaseContagion::seed,
                R"pbdoc(
            Seed the RNG.

            Args:
               seed: seed for the RNG.
            )pbdoc", py::arg("seed"))

        .def("evolve", &BaseContagion::evolve,
                R"pbdoc(
            Let the system evolve over a period of time.

            Args:
               period: Time period of the evolution.
               decorrelation_time (optional): Time period for decorrelation.
               measure: Bool, if true measure the system after decorrelation.
               quasistationary: Bool, if true quasistationary state.
            )pbdoc", py::arg("period"), py::arg("decorrelation_time")=1,
                py::arg("measure")=false, py::arg("quasistationary")=false)


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
            )pbdoc")

        .def("measure_infectious_set",
        &MeasurableContagionProcess::measure_infectious_set,
            R"pbdoc(
            Add the measure of infectious set to the vector
            of measures to be performed during simulations.
            )pbdoc")

        .def("measure_time",
        &MeasurableContagionProcess::measure_time,
            R"pbdoc(
            Add the measure of time to the vector
            of measures to be performed during simulations.
            )pbdoc");


    /* =================================
     * Class deriving from BaseContagion
     * =================================*/

    py::class_<GroupSIS, BaseContagion>(m, "GroupSIS")

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

        .def("get_lifetime", &GroupSIS::get_lifetime, R"pbdoc(
            Returns the lifetime for the current state.
            )pbdoc");

    py::class_<ContinuousSIS, BaseContagion>(m, "ContinuousSIS")

        .def(py::init<EdgeList&, double,
                const vector<vector<double>>&,
                const vector<double>&>(), R"pbdoc(
            Default constructor of the class ContinuousSIS.

            Args:
               edge_list: Edge list for the network structure.
               recovery_rate: Double for the recovery rate
               infection_rate: Matrix of infection rate per node in groups
               group_transmission_rate: Vector of transmission rate per group
            )pbdoc", py::arg("edge_list"),
                py::arg("recovery_rate"),
                py::arg("infection_rate"),
                py::arg("group_transmission_rate"))

        .def("get_lifetime", &ContinuousSIS::get_lifetime, R"pbdoc(
            Returns the lifetime for the current state.
            )pbdoc");


    py::class_<ContinuousSIR, BaseContagion>(m, "ContinuousSIR")

        .def(py::init<EdgeList&, double,
                const vector<vector<double>>&,
                const vector<double>&>(), R"pbdoc(
            Default constructor of the class ContinuousSIR.

            Args:
               edge_list: Edge list for the network structure.
               recovery_rate: Double for the recovery rate
               infection_rate: Matrix of infection rate per node in groups
               group_transmission_rate: Vector of transmission rate per group
            )pbdoc", py::arg("edge_list"),
                py::arg("recovery_rate"),
                py::arg("infection_rate"),
                py::arg("group_transmission_rate"))

        .def("get_lifetime", &ContinuousSIR::get_lifetime, R"pbdoc(
            Returns the lifetime for the current state.
            )pbdoc");


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

    py::class_<DiscreteSIS, BaseContagion>(m, "DiscreteSIS")

        .def(py::init<EdgeList, double,
                std::vector<std::vector<double>>>(), R"pbdoc(
            Default constructor of the class GroupSIS.

            Args:
               edge_list: Edge list for the network structure.
               recovery_probability: Double for the recovery probability
               infection_probability: vector of vector for the infection
                                      probability for different group size
                                      and number of infected
            )pbdoc", py::arg("edge_list"),
                py::arg("recovery_probability"),
                py::arg("infection_probability"))

        .def("get_lifetime", &DiscreteSIS::get_lifetime, R"pbdoc(
            Returns the lifetime for the current state.
            )pbdoc");


    py::class_<HeterogeneousExposure, BaseContagion>(m, "HeterogeneousExposure")

        .def(py::init<EdgeList, double, double, double,
                double, double>(), R"pbdoc(
            Default constructor of the class GroupSIS.

            Args:
               edge_list: Edge list for the network structure.
               recovery_probability: Double for the recovery probability
               alpha:  Double for the exponent of the participation time distribution
               T: Double for the temporal window
               beta: Double for the rate of dose accumulation
               K: Double for the dose threshold
            )pbdoc", py::arg("edge_list"),
                py::arg("recovery_probability"),
                py::arg("alpha"),
                py::arg("T"),
                py::arg("beta"),
                py::arg("K"))

        .def("get_lifetime", &HeterogeneousExposure::get_lifetime, R"pbdoc(
            Returns the lifetime for the current state.
            )pbdoc");



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

    py::class_<InfectiousSet,shared_ptr<InfectiousSet>>(m,
                "InfectiousSet")

        .def(py::init<>(), R"pbdoc(
            Default constructor of the class InfectiousSet.

            )pbdoc")

        .def("get_name", &InfectiousSet::get_name, R"pbdoc(
            Returns the name of the measure.
            )pbdoc")

        .def("get_result", &InfectiousSet::get_result, R"pbdoc(
            Returns the result associated to the measure.
            )pbdoc");

    py::class_<Time,shared_ptr<Time>>(m,
                "Time")

        .def(py::init<>(), R"pbdoc(
            Default constructor of the class Time.

            )pbdoc")

        .def("get_name", &Time::get_name, R"pbdoc(
            Returns the name of the measure.
            )pbdoc")

        .def("get_result", &Time::get_result, R"pbdoc(
            Returns the result associated to the measure.
            )pbdoc");

}
