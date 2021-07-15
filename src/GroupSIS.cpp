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
    BaseContagion(edge_list),
    recovery_rate_(recovery_rate), infection_rate_(infection_rate),
    event_set_(rate_bounds.first,rate_bounds.second)
{
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


//advance the process to the next step by performing infection/recovery
//it is assumed that the lifetime is finite
inline void GroupSIS::next_event()
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
void GroupSIS::clear()
{
    BaseContagion::clear();
    event_set_.clear(); //to avoid numerical error accumulation
}



}//end of namespace schon
