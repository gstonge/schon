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

#include "ContinuousSIR.hpp"
#include <optional>
#include <utility>
#include <algorithm>
#include <iostream>
#include <exception>
#include <limits>

using namespace std;

namespace schon
{//start of namespace schon

//constructor of the class
ContinuousSIR::ContinuousSIR(const EdgeList& edge_list, double recovery_rate,
        const vector<vector<double>>& infection_rate,
        const vector<double>& group_transmission_rate):
    BaseContagion(edge_list),
    recovery_rate_(recovery_rate),
    infection_rate_(infection_rate),
    group_transmission_rate_(group_transmission_rate),
    event_set_(1.,1.)
{
    //determine min/max rate upper and lower bounds
    double min_transmission = std::numeric_limits<double>::infinity();
    double max_transmission = 0;
    for (double rate : group_transmission_rate_)
    {
        if (rate > 0)
        {
            if (rate < min_transmission)
            {
                min_transmission = rate; //min non zero transmission
            }
            if (rate > max_transmission)
            {
                max_transmission = rate;
            }
        }
    }
    double min = recovery_rate_;
    double max = recovery_rate_;
    for (int n = 2; n < infection_rate_.size(); n++)
    {
        for (int i = 0; i <= n; i++)
        {
            double rate = (n-i)*infection_rate_[n][i];
            if (rate > 0)
            {
                if (min_transmission*rate < min)
                {
                    min = min_transmission*rate;
                }
                if (max_transmission*rate > max)
                {
                    max = max_transmission*rate;
                }
            }
        }
    }
    event_set_ = sset::SamplableSet<Event>(min,max); //set true bounds

}

//update the event group rate
inline void ContinuousSIR::update_group_rate(Group group, Node node,
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
inline void ContinuousSIR::infect(Node node)
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
inline void ContinuousSIR::recover(Node node)
{
    if (node_state_vector_[node] == I)
    {
        node_state_vector_[node] = R;
        infected_node_set_.erase(node);
        for (Group group : network_.adjacent_groups(node))
        {
            update_group_rate(group,node,I,R);
        }
        //erase the recovery event for the node
        event_set_.erase(make_tuple(NODE,RECOVERY,node));
    }
    else
    {
        throw runtime_error("Recovery attempt: the node is not infected");
    }
}


//advance the process to the next step by performing infection/recovery
//it is assumed that the lifetime is finite
inline void ContinuousSIR::next_event()
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



//clear the state; as if all node became susceptible at this time
//clear all measures as well
//overload BaseContagion
void ContinuousSIR::clear()
{
    BaseContagion::clear();
    event_set_.clear(); //to avoid numerical error accumulation
}



}//end of namespace schon
