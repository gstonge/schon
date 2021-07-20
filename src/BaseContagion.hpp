/*
 * MIT License
 *
 * Copyright (c) 2020 Guillaume St-Onge
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

#ifndef BASECONTAGION_HPP_
#define BASECONTAGION_HPP_

#include "MeasurableContagionProcess.hpp"
#include "SamplableSet/SamplableSet.hpp"
#include <iostream>

namespace schon
{//start of namespace schon


//abstract class with more functionality to avoid overlapp between classes
class BaseContagion : public MeasurableContagionProcess
{
public:
    //Constructor
    BaseContagion(const EdgeList& edge_list);

    //Accessors
    std::size_t size() const
        {return network_.size();}
    const std::vector<NodeState>& get_node_state_vector() const
        {return node_state_vector_;}
    const std::unordered_set<Node>& get_infected_node_set() const
        {return infected_node_set_;}
    const BipartiteNetwork& get_network() const
        {return network_;}
    double get_current_time() const
        {return current_time_;}
    std::size_t get_number_of_infected_nodes() const
        {return infected_node_set_.size();}
    double get_lifetime() const //dummy definition
        {return 1;}

    //Mutators
    void seed(unsigned int seed)
        {gen_.seed(seed);}
    void infect_fraction(double fraction);
    void infect_node_set(const std::unordered_set<Node>& node_set);

    void clear();
    void reset();
    void initialize_history(std::size_t number_of_states = 100);

    void evolve(double period, double decorrelation_time=1, bool measure=false,
            bool quasistationary=false);


    //it seems abstract base class definition in the binding does not work

protected:
    //Members
    BipartiteNetwork network_;
    std::vector<NodeState> node_state_vector_;
    std::vector<GroupState> group_state_vector_;
    std::vector<GroupStatePosition> group_state_position_vector_;
    std::unordered_set<Node> infected_node_set_;
    std::vector<std::unordered_set<Node>> history_vector_;

    double current_time_;
    double last_event_time_;
    double time_since_last_measure_;
    sset::RNGType &gen_;
    mutable std::uniform_real_distribution<double> random_01_;

    //utility functions
    Node random_node(Group group, NodeState node_state) const;
    void store_configuration();
    void get_configuration_from_history();

    void infect(Node node) {}; //dummy definition
    void recover(Node node) {}; //dummy definition
    void next_event() {}; //dummy definition


};

}//end of namespace schon

#endif /* BASECONTAGION_HPP_ */
