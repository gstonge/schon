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

#ifndef HETEROGENEOUSEXPOSURE_HPP_
#define HETEROGENEOUSEXPOSURE_HPP_

#include "BaseContagion.hpp"

namespace schon
{//start of namespace schon


//class to simulate SIS process on networks
class HeterogeneousExposure : public BaseContagion
{
public:
    //Constructor
    HeterogeneousExposure(const EdgeList& edge_list, double recovery_probability,
            double alpha, double T, double beta, double K);

    //Accessors
    double get_lifetime() const
        {return infected_node_set_.size() == 0 ?
            std::numeric_limits<double>::infinity() : 1.;}

    //Mutators
    void clear();

protected:
    //Members
    double recovery_probability_;
    double recovery_propensity_;
    sset::SamplableSet<Node> recovery_event_set_;
    std::poisson_distribution<int> poisson_dist_;
    double alpha_;
    double T_;
    double beta_;
    double K_;

    //utility functions

    inline double get_dose(double tau, double rho)
    {
        double r = random_01_(gen_);
        return -beta_*tau*rho*log(1-r);
    }
    inline double get_participation_time()
    {
        double r = random_01_(gen_);
        return pow(1/(1-r*(1-pow(T_,-alpha_))), 1./alpha_);
    }

    inline void update_group_state(Group group, Node node,
            NodeState previous_state, NodeState new_state);

    inline void infect(Node node);
    inline void recover(Node node);
    inline void next_event();
};

}//end of namespace schon

#endif /* HETEROGENEOUSEXPOSURE_HPP_ */

