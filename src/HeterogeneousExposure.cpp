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

#include "HeterogeneousExposure.hpp"
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
HeterogeneousExposure::HeterogeneousExposure(const EdgeList& edge_list, double recovery_probability,
        double alpha, double T, double beta, double K):
    BaseContagion(edge_list),
    recovery_probability_(recovery_probability),
    recovery_propensity_(-log(1-recovery_probability)),
    recovery_event_set_(1.,1.),
    poisson_dist_(1.),
    alpha_(alpha),
    T_(T),
    beta_(beta),
    K_(K)
{
}

//update the group state
inline void HeterogeneousExposure::update_group_state(Group group, Node node,
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
}

//infect a node
inline void HeterogeneousExposure::infect(Node node)
{
    if (node_state_vector_[node] == S)
    {
        node_state_vector_[node] = I;
        infected_node_set_.insert(node);
        for (Group group : network_.adjacent_groups(node))
        {
            update_group_state(group,node,S,I);
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
inline void HeterogeneousExposure::recover(Node node)
{
    if (node_state_vector_[node] == I)
    {
        node_state_vector_[node] = S;
        infected_node_set_.erase(node);
        for (Group group : network_.adjacent_groups(node))
        {
            update_group_state(group,node,I,S);
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
inline void HeterogeneousExposure::next_event()
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
    //get the infections
    double tau, kappa, rho;
    double n,i;
    unordered_set<Node> new_infected;
    for (Group group = 0; group < network_.number_of_groups(); group++)
    {
        GroupState & group_state = group_state_vector_[group];
        n = network_.group_size(group);
        i = group_state[I].size();
        rho = i/(n-1);
        //for all susceptible, check for infections
        for (Node node : group_state[S])
        {
            tau = get_participation_time();
            kappa = get_dose(tau,rho);
            if (kappa > K_)
            {
                new_infected.insert(node);
            }
        }
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
void HeterogeneousExposure::clear()
{
    BaseContagion::clear();
    recovery_event_set_.clear(); //to avoid numerical error accumulation
}



}//end of namespace schon
