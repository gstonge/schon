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

#include "DiscreteSIS.hpp"
#include <optional>
#include <utility>
#include <iostream>
#include <exception>
#include <cmath>
#include <limits>

using namespace std;

namespace schon
{//start of namespace schon

//constructor of the class
DiscreteSIS::DiscreteSIS(const EdgeList& edge_list, double recovery_probability,
        const std::vector<std::vector<double>>& infection_probability):
    BaseContagion(edge_list),
    recovery_probability_(recovery_probability),
    recovery_propensity_(-log(1-recovery_probability)),
    infection_probability_(infection_probability),
    infection_propensity_(infection_probability.size(),vector<double>()),
    infection_event_set_(1.,1.),
    recovery_event_set_(1.,1.),
    poisson_dist_(1.)
{
    //calculate Poisson rate equivalent for each probability
    double min = std::numeric_limits<double>::infinity();
    double max = 0;
    double propensity = 0;
    for (int i = 0; i < infection_probability.size(); i++)
    {
        for (auto prob : infection_probability[i])
        {
            propensity = -log(1-prob);
            infection_propensity_[i].push_back(propensity);
            if (propensity > max)
            {
                max = propensity;
            }
            if (propensity < min and propensity > 0)
            {
                min = propensity;
            }
        }
    }
    max *= network_.max_group_size(); //upper bound
    infection_event_set_ = sset::SamplableSet<Group>(min,max); //set true bounds
}

//update the group state and the infection propensity
inline void DiscreteSIS::update_infection_propensity(Group group, Node node,
        NodeState previous_state, NodeState new_state)
{
    GroupState & group_state = group_state_vector_[group];
    size_t position = group_state_position_vector_[group][node];
    //if node not in the back, put node in the back position
    swap(group_state[previous_state][position],
            group_state[previous_state].back());
    //also update the position of the node in the back
    Node back_node = group_state[previous_state][position];
    group_state_position_vector_[group][back_node] = position;
    //put node in new group state
    group_state[previous_state].pop_back();
    group_state_position_vector_[group][node] = group_state[new_state].size();
    group_state[new_state].push_back(node);
    //update event set with new propensity
    double new_propensity = get_infection_propensity(group);
    if (new_propensity > 0)
    {
        infection_event_set_.set_weight(group,new_propensity);
    }
    else
    {
        infection_event_set_.erase(group);
    }
}

//infect a node
inline void DiscreteSIS::infect(Node node)
{
    if (node_state_vector_[node] == S)
    {
        node_state_vector_[node] = I;
        infected_node_set_.insert(node);
        for (Group group : network_.adjacent_groups(node))
        {
            update_infection_propensity(group,node,S,I);
        }
        //create a recovery event for the node
        recovery_event_set_.insert(node, 1.);
    }
    else
    {
        throw runtime_error("Infection attempt: the node is not susceptible");
    }
}

//recover a node
inline void DiscreteSIS::recover(Node node)
{
    if (node_state_vector_[node] == I)
    {
        node_state_vector_[node] = S;
        infected_node_set_.erase(node);
        for (Group group : network_.adjacent_groups(node))
        {
            update_infection_propensity(group,node,I,S);
        }
        //erase the recovery event for the node
        recovery_event_set_.erase(node);
    }
    else
    {
        throw runtime_error("Recovery attempt: the node is not infected");
    }
}


//advance the process to the next step by performing infection/recovery
inline void DiscreteSIS::next_event()
{
    current_time_ = last_event_time_ + get_lifetime();
    //get the number of recoveries and assign them
    poisson_dist_ = poisson_distribution<int>(
            recovery_propensity_*recovery_event_set_.size());
    int nb_rec = poisson_dist_(gen_);
    unordered_set<Node> new_susceptible; //use a set to discard repetition
    for (int i = 0; i < nb_rec; i++)
    {
        pair<Node,double> node_weight_pair =
            (recovery_event_set_.sample()).value();
        new_susceptible.insert(node_weight_pair.first);
    }
    //get the number of infections and assign them
    poisson_dist_ = poisson_distribution<int>(
            infection_event_set_.total_weight());
    int nb_inf = poisson_dist_(gen_);
    unordered_set<Node> new_infected;
    for (int i = 0; i < nb_inf; i++)
    {
        pair<Group,double> group_weight_pair =
            (infection_event_set_.sample()).value();
        Group group = group_weight_pair.first;
        //choose uniformly among susceptible
        new_infected.insert(random_node(group, S));
    }
    //perform recovery and infections
    for (Node node : new_susceptible)
    {
        recover(node);
    }
    for (Node node : new_infected)
    {
        infect(node);
    }
    last_event_time_ = current_time_;
}


//clear the state; as if all node became susceptible at this time
//clear all measures as well
//overload BaseContagion
void DiscreteSIS::clear()
{
    BaseContagion::clear();
    infection_event_set_.clear(); //to avoid numerical error accumulation
    recovery_event_set_.clear(); //to avoid numerical error accumulation
}



}//end of namespace schon
