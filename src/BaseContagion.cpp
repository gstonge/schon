/*
 * MIT License
 *
 * Copyright (c) 2021 Guillaume St-Onge
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

#include "BaseContagion.hpp"
#include <optional>
#include <utility>
#include <iostream>
#include <exception>

using namespace std;

namespace schon
{//start of namespace schon

//constructor of the class
BaseContagion::BaseContagion(const EdgeList& edge_list):
    network_(edge_list),
    node_state_vector_(network_.size(), S),
    group_state_vector_(network_.number_of_groups()),
    group_state_position_vector_(network_.number_of_groups()),
    infected_node_set_(),
    history_vector_(),
    current_time_(0),
    last_event_time_(0),
    time_since_last_measure_(0),
    gen_(sset::BaseSamplableSet::gen_),
    random_01_()
{
    //initialize group state vector and position
    for (Group group : network_.groups())
    {
        for (unsigned int state = S; state < STATECOUNT; state++)
        {
            group_state_vector_[group].push_back(vector<Node>());
        }
        for (Node node : network_.group_members(group))
        {
            group_state_position_vector_[group][node] =
                group_state_vector_[group][0].size();
            group_state_vector_[group][0].push_back(node); //nodes are S
        }
    }
}

//get a random node of the particular state in the group
Node BaseContagion::random_node(Group group, NodeState node_state) const
{
    const GroupState& group_state = group_state_vector_[group];
    unsigned int index = floor(random_01_(gen_)*group_state[node_state].size());
    return group_state[node_state][index];
}

//infect a fraction of the nodes
void BaseContagion::infect_fraction(double fraction)
{
    unsigned int number_of_infection = floor(network_.size()*fraction);
    Node node;
    unsigned int count = 0;
    while (count < number_of_infection)
    {
        node = floor(random_01_(gen_)*network_.size());
        if (node_state_vector_[node] == S)
        {
            infect(node);
            count += 1;
        }
    }
}

//infect a certain set of of nodes
void BaseContagion::infect_node_set(const std::unordered_set<Node>& node_set)
{
    for (Node node : node_set)
    {
        if (node_state_vector_[node] == S)
        {
            infect(node);
        }
    }
}


//clear the state; as if all node became susceptible at this time
void BaseContagion::clear()
{
    //recover nodes
    for (Node node : infected_node_set_)
    {
        recover(node);
    }
}

//clear and reset the process to initial state at time 0 (and clear history)
//clear all measures as well
void BaseContagion::reset()
{
    clear();
    //clear measures
    for(size_t i = 0; i < measure_vector_.size(); i++)
    {
        measure_vector_[i] -> clear();
    }
    history_vector_.clear();
    current_time_ = 0;
    last_event_time_ = 0;
    time_since_last_measure_ = 0;
}

void BaseContagion::initialize_history(std::size_t number_of_states)
{
    if (history_vector_.size() > 0)
    {
        history_vector_.clear();
    }
    //must be non trivial
    history_vector_ = std::vector<std::unordered_set<Node>>(number_of_states,
            infected_node_set_);
}

void BaseContagion::store_configuration()
{
    //delete randomly one configuration before
    size_t index = floor(random_01_(gen_)*history_vector_.size());
    swap(history_vector_[index], history_vector_.back());
    history_vector_.pop_back();
    history_vector_.push_back(infected_node_set_);
}

void BaseContagion::get_configuration_from_history()
{
    clear();
    size_t index = floor(random_01_(gen_)*history_vector_.size());
    const std::unordered_set<Node>& infected_node_set = history_vector_[index];
    for (Node node : infected_node_set)
    {
        infect(node);
    }
}

//perform the evolution of the process over a period of time and perform
//measures after each decorrelation time if needed
void BaseContagion::evolve(double period, double decorrelation_time, bool measure,
        bool quasistationary)
{
    if (quasistationary and history_vector_.size() == 0)
    {
        initialize_history();
    }
    double initial_time = current_time_;
    while(last_event_time_ + get_lifetime() - initial_time < period)
    {
        time_since_last_measure_ += last_event_time_
            + get_lifetime() - current_time_; //after the coming event
        if (time_since_last_measure_ > decorrelation_time)
        {
            time_since_last_measure_ -= decorrelation_time;
            if (measure)
            {
                for(size_t i = 0; i < measure_vector_.size(); i++)
                {
                    measure_vector_[i] -> measure(this);
                }
            }
            if (quasistationary)
            {
                store_configuration();
            }
        }
        next_event();
        if (isinf(get_lifetime()) and quasistationary)
        {
            get_configuration_from_history();
        }
    }
    time_since_last_measure_ += period - (last_event_time_ - initial_time);
    //if we need to perform a last measure
    if (time_since_last_measure_ > decorrelation_time)
    {
        time_since_last_measure_ -= decorrelation_time;
        if (measure)
        {
            for(size_t i = 0; i < measure_vector_.size(); i++)
            {
                measure_vector_[i] -> measure(this);
            }
        }
        if (quasistationary)
        {
            store_configuration();
        }
    }
    current_time_ = initial_time + period;
}



}//end of namespace schon
