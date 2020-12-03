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

#include "ComplexSIS.hpp"
#include <optional>
#include <utility>

using namespace std;

namespace schon
{//start of namespace schon

//constructor of the class
ComplexSIS::ComplexSIS(const EdgeList& edge_list,
        const function<double(size_t,size_t)>& recovery_rate,
        const function<double(size_t,size_t)>& infection_rate,
        const pair<double>& group_rate_bounds):
    current_time_(0), last_event_time_(0), time_since_last_measure_(0),
    gen_(sset::BaseSamplableSet::gen_), random_01_(),
    recovery_rate_(recovery_rate), infection_rate_(infection_rate),
    network_(edge_list), node_state_vector_(network_.size(), S),
    group_state_vector_(network_.number_of_groups()),
    group_state_position_vector_(network_.number_of_groups()),
    infected_node_set(),
    event_set(group_rate_bounds.first,group_rate_bounds.second)
{
    //initialize group state vector and position
    for (Group group : network_.groups())
    {
        for (NodeState state = 0; state < STATECOUNT; state++)
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

//infect a node
inline void ComplexSIS::infect(Node node)
{
    NodeState previous_state = state_vector_[node];
    state_vector_[node] = I;
    infected_node_set_.insert(node);
    for (Group group : network_.adjacent_groups(node))
    {
        //update group state and group state position of node
        GroupState & group_state = group_state_vector_[group];
        size_t position = group_state_position_vector_[group][node];
        swap(group_state[previous_state][position],group_state.back());
        group_state[previous_state].pop_back();
        group_state_position_vector_[group][node] = group_state[I].size();
        group_state[I].push_back(node);
        //update event set
        double new_rate = get_recovery_rate(group) +
            get_infection_rate(group);
        if (new_rate > 0)
        {
            event_set.set_weight(group,new_rate);
        }
        else
        {
            event_set.erase(group);
        }
    }
}

//recover a node
inline void ComplexSIS::recover(Node node)
{
    NodeState previous_state = state_vector_[node];
    state_vector_[node] = S;
    infected_node_set_.erase(node);
    for (Group group : network_.adjacent_groups(node))
    {
        //update group state and group state position of node
        GroupState & group_state = group_state_vector_[group];
        size_t position = group_state_position_vector_[group][node];
        swap(group_state[previous_state][position],group_state.back());
        group_state[previous_state].pop_back();
        group_state_position_vector_[group][node] = group_state[S].size();
        group_state[S].push_back(node);
        //update event set
        double new_rate = get_recovery_rate(group) +
            get_infection_rate(group);
        if (new_rate > 0)
        {
            event_set.set_weight(group,new_rate);
        }
        else
        {
            event_set.erase(group);
        }
    }
}

//infect a fraction of the nodes
void ComplexSIS::infect_fraction(double fraction)
{
	unsigned int number_of_infection = floor(network_.size()*fraction);
	Node node;
    unsigned int count = 0;
	while (count < number_of_infection)
	{
		node = floor(random_01_(gen_)*network_.size());
		if (state_vector_[node] == S)
		{
			infect(node);
            count += 1;
		}
	}
}

//get a random node of the particular state in the group
Node ComplexSIS::random_node(Group group, NodeState node_state) const
{
    GroupState& group_state = group_state_vector_[group];
    unsigned int index = floor(random_01_(gen_)*group_state[node_state].size());
    return group_state[node_state][index];
}


//advance the process to the next step by performing infection/recovery
//it is assumed that the lifetime is finite
void ComplexSIS::next_event()
{
    current_time_ = last_event_time_ + get_lifetime();
    //select a group proportionally to its weight
    pair<Group, double> group_weight_pair = (event_set.sample()).value();
    const Group& group= group_weight_pair.first;
    if(random_01_(gen_) < get_recovery_rate(group)/group_weight_pair.second)
    {
        //recovery event
        Node node = random_node(group, I);
        recover(node);
    }
    else
    {
        //infection event
        Node node = random_node(group, S);
        infect(node);
    }
    last_event_time_ = current_time_;
}

//perform the evolution of the process over a period of time
void ComplexSIS::evolve(double period)
{
    double initial_time = current_time_;
    while(last_event_time_ + get_lifetime() - initial_time < period)
    {
        next_event();
    }
    current_time_ = initial_time + period;
}

//perform the evolution of the process over a period of time and perform
//measures after each decorrelation time
void ComplexSIS::evolve_and_measure(double period, double decorrelation_time)
{
    double initial_time = current_time_;
    while(last_event_time_ + get_lifetime() - initial_time < period)
    {
        time_since_last_measure_ += last_event_time_
            + get_lifetime() - current_time_; //after the coming event
        if (time_since_last_measure_ > decorrelation_time)
        {
            time_since_last_measure_ -= decorrelation_time;
            for(size_t i = 0; i < measure_vector_.size(); i++)
            {
                measure_vector_[i] -> measure(this);
            }
        }
        next_event();
    }
    time_since_last_measure_ += period - (last_event_time_ - initial_time);
    //if we need to perform a last measure
    if (time_since_last_measure_ > decorrelation_time)
    {
        for(size_t i = 0; i < measure_vector_.size(); i++)
        {
            measure_vector_[i] -> measure(this);
        }
    }
    current_time_ = initial_time + period;
}

//clear the state; as if all node became susceptible at this time
void ComplexSIS::clear()
{
    for (Node node : infected_node_set_)
    {
        recover(node)
    }
    event_set.clear();
}

//clear and reset the process to initial state at time 0
void ComplexSIS::reset()
{
    clear();
    current_time_ = 0;
    last_event_time_ = 0;
    time_since_last_measure_ = 0;
}

////set the state of the process using an SIS configuration--keep current time
//void ComplexSIS::set_configuration(
        //const SIS_configuration& configuration)
//{
    //clear();
    ////time shift between current configuration and the one in argument
    //double time_shift = current_time_ - configuration.get_current_time();
    //const unordered_set<Node>& infected_node_set =
        //configuration.get_infected_node_set();
    //const unordered_map<Node,double>& time_of_infection_map =
        //configuration.get_time_of_infection_map();
    //for (auto iter = infected_node_set.begin();
            //iter != infected_node_set.end(); iter ++)
    //{
        //delayed_infect(*iter, time_shift + time_of_infection_map.at(*iter));
    //}
//}

}//end of namespace schon
