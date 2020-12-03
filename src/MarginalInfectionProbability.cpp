/*
 * MIT License
 *
 * Copyright (c) 2018 Guillaume St-Onge
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

#include "MarginalInfectionProbability.hpp"

using namespace std;

namespace spreadnet
{//start of namespace spreadnet

//constructor
MarginalInfectionProbability::MarginalInfectionProbability(
        size_t network_size): weight_vector_(network_size, 0.),
    count_(0), name_("marginal_infection_probability")
{
}

//return the marginal probability of infection for each node
vector<double> MarginalInfectionProbability::get_result() const
{
    vector<double> marginal_vector(weight_vector_);
    if (count_ > 0)
    {
        for (auto iter = marginal_vector.begin();
                iter != marginal_vector.end(); iter++)
        {
            *iter /= count_;
        }
    }
    return marginal_vector;
}

//perform a measure on the spreading process
void MarginalInfectionProbability::measure(
        SpreadingProcess const * const ptr)
{
    const unordered_set<NodeLabel>& infected_node_set =
        ptr->get_infected_node_set();
    //iterate on infected nodes
    for (auto iter = infected_node_set.begin();
            iter != infected_node_set.end(); iter++)
    {
        weight_vector_[*iter] += 1;
    }
    count_ += 1;
}

}//end of namespace spreadnet
