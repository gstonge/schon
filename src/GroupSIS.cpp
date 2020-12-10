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

#include "GroupSIS.hpp"
#include <optional>
#include <utility>
#include <iostream>
#include <exception>

using namespace std;

namespace schon
{//start of namespace schon

//constructor of the class
GroupSIS::GroupSIS(const EdgeList& edge_list, double recovery_rate,
        const function<double(size_t,size_t)>& infection_rate,
        const pair<double,double>& rate_bounds):
    current_time_(0), last_event_time_(0), time_since_last_measure_(0),
    gen_(sset::BaseSamplableSet::gen_), random_01_(),
    recovery_rate_(recovery_rate), infection_rate_(infection_rate),
    network_(edge_list), node_state_vector_(network_.size(), S),
    group_state_vector_(network_.number_of_groups()),
    group_state_position_vector_(network_.number_of_groups()),
    infected_node_set_(),
    event_set_(rate_bounds.first,rate_bounds.second),
    history_vector_()
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

//update the event group rate
inline void GroupSIS::update_group_rate(Group group, Node node,
        NodeState previous_state, NodeState new_state)
{
    GroupState & group_state = group_state_vector_[group];
    size_t position = group_state_position_vector_[group][node];
    if (position < group_state[previous_state].size())
    {
        //node not in the back, put node in the back position
        swap(group_state[previous_state][position],
                group_state[previous_state].back());
        //also update the position of the node in the back
        Node back_node = group_state[previous_state][position];
        group_state_position_vector_[group][back_node] = position;
    }
    //put node in new group state
    group_state[previous_state].pop_back();
    group_state_position_vector_[group][node] = group_state[new_state].size();
    group_state[new_state].push_back(node);
    //update event set with new rate
    double new_rate = get_infection_rate(group);
    if (new_rate > 0)
    {
        event_set_.set_weight(make_tuple(GROUP,INFECTION,group),new_rate);
    }
    else
    {
        event_set_.erase(make_tuple(GROUP,INFECTION,group));
    }
}

//infect a node
inline void GroupSIS::infect(Node node)
{
    if (node_state_vector_[node] == S)
    {
        node_state_vector_[node] = I;
        infected_node_set_.insert(node);
        for (Group group : network_.adjacent_groups(node))
        {
            update_group_rate(group,node,S,I);
        }
        //create a recovery event for the node
        event_set_.insert(make_tuple(NODE,RECOVERY,node), recovery_rate_);
    }
    else
    {
        throw runtime_error("Infection attempt: the node is not susceptible");
    }
}

//recover a node
inline void GroupSIS::recover(Node node)
{
    if (node_state_vector_[node] == I)
    {
        node_state_vector_[node] = S;
        infected_node_set_.erase(node);
        for (Group group : network_.adjacent_groups(node))
        {
            update_group_rate(group,node,I,S);
        }
        //erase the recovery event for the node
        event_set_.erase(make_tuple(NODE,RECOVERY,node));
    }
    else
    {
        throw runtime_error("Recovery attempt: the node is not infected");
    }
}

//infect a fraction of the nodes
void GroupSIS::infect_fraction(double fraction)
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

//get a random node of the particular state in the group
inline Node GroupSIS::random_node(Group group, NodeState node_state) const
{
    const GroupState& group_state = group_state_vector_[group];
    unsigned int index = floor(random_01_(gen_)*group_state[node_state].size());
    return group_state[node_state][index];
}


//advance the process to the next step by performing infection/recovery
//it is assumed that the lifetime is finite
void GroupSIS::next_event()
{
    current_time_ = last_event_time_ + get_lifetime();
    //select a group proportionally to its weight
    pair<Event, double> event_weight_pair = (event_set_.sample()).value();
    const Event& event = event_weight_pair.first;
    if (get<0>(event) == NODE and get<1>(event) == RECOVERY)
    {
        //node-based recovery event
        Node node = get<2>(event);
        recover(node);
    }
    else if (get<0>(event) == GROUP and get<1>(event) == INFECTION)
    {
        //Groub-based infection event
        Group group = get<2>(event);
        Node node = random_node(group, S);
        infect(node);
    }
    else
    {
        throw runtime_error("Unallowed type of event");
    }
    last_event_time_ = current_time_;
}


//perform the evolution of the process over a period of time and perform
//measures after each decorrelation time if needed
void GroupSIS::evolve(double period, double decorrelation_time, bool measure,
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

//clear the state; as if all node became susceptible at this time
void GroupSIS::clear()
{
    for (Node node : infected_node_set_)
    {
        recover(node);
    }
    event_set_.clear();
}

//clear and reset the process to initial state at time 0 (and clear history)
void GroupSIS::reset()
{
    clear();
    history_vector_.clear();
    current_time_ = 0;
    last_event_time_ = 0;
    time_since_last_measure_ = 0;
}

void GroupSIS::initialize_history(std::size_t number_of_states)
{
    if (history_vector_.size() > 0)
    {
        history_vector_.clear();
    }
    //must be non trivial
    history_vector_ = std::vector<std::unordered_set<Node>>(number_of_states,
            infected_node_set_);
}

void GroupSIS::store_configuration()
{
    //delete randomly one configuration before
    size_t index = floor(random_01_(gen_)*history_vector_.size());
    swap(history_vector_[index], history_vector_.back());
    history_vector_.pop_back();
    history_vector_.push_back(infected_node_set_);
}

void GroupSIS::get_configuration_from_history()
{
    clear();
    size_t index = floor(random_01_(gen_)*history_vector_.size());
    const std::unordered_set<Node>& infected_node_set = history_vector_[index];
    for (Node node : infected_node_set)
    {
        infect(node);
    }
}


}//end of namespace schon
