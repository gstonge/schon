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

#ifndef COMPLEXSIS_HPP_
#define COMPLEXSIS_HPP_

#include <SamplableSet.hpp>
#include "MeasurableContagionProcess.hpp"
#include <functional>

namespace schon
{//start of namespace schon


//class to simulate SIS process on networks
class ComplexSIS : public MeasurableContagionProcess
{
public:
    //Constructor
    ComplexSIS(const EdgeList& edge_list,
            const std::function<double(std::size_t,std::size_t)>& recovery_rate,
            const std::function<double(std::size_t,std::size_t)>& infection_rate,
            const std::pair<double,double>& group_rate_bounds);

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
    double get_lifetime() const
        {return event_set_.size() == 0 ?
            std::numeric_limits<double>::infinity() :
            1/event_set_.total_weight();}
    std::size_t get_number_of_infected_nodes() const
        {return infected_node_set_.size();}

    double get_recovery_rate(Group group) const
        {return recovery_rate_(network_.group_size(group),
                group_state_vector_[group][I].size());}
    double get_infection_rate(Group group) const
        {return infection_rate_(network_.group_size(group),
                group_state_vector_[group][I].size());}

    //Mutators
    void seed(unsigned int seed)
        {gen_.seed(seed);}
    void infect_fraction(double fraction);
    void next_event();
    void evolve(double period);
    void evolve_and_measure(double period, double decorrelation_time);
    void clear();
    void reset();

protected:
    double current_time_;
    double last_event_time_;
    double time_since_last_measure_;
    sset::RNGType &gen_;
    mutable std::uniform_real_distribution<double> random_01_;
    std::function<double(std::size_t,std::size_t)> recovery_rate_;
    std::function<double(std::size_t,std::size_t)> infection_rate_;

private:
    //Members
    BipartiteNetwork network_;
    std::vector<NodeState> node_state_vector_;
    std::vector<GroupState> group_state_vector_;
    std::vector<GroupStatePosition> group_state_position_vector_;
    std::unordered_set<Node> infected_node_set_;
    sset::SamplableSet<Group> event_set_;

    //utility functions
    inline void update_group_rate(Group group, Node node,
            NodeState previous_state, NodeState new_state);
    inline void infect(Node node);
    inline void recover(Node node);
    inline Node random_node(Group group, NodeState node_state) const;
};

}//end of namespace schon

#endif /* COMPLEXSIS_HPP_ */

