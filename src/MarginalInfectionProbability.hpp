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

#ifndef MARGINALINFECTIONPROBABILITY_HPP_
#define MARGINALINFECTIONPROBABILITY_HPP_

#include "Measure.hpp"

namespace schon
{//start of namespace schon

class MarginalInfectionProbability : public Measure
{
public:
    //Constructor
    MarginalInfectionProbability(std::size_t network_size);

    //Acessors
    std::vector<double> get_result() const;
    const std::string& get_name() const
        {return name_;}

    //Mutators
    void measure(ContagionProcess const * const ptr);
    void clear()
        {count_ = 0; weight_vector_ = std::vector<double>(network_size, 0.);}

private:
    //Members
    const std::string name_;
    int count_;
    std::vector<double> weight_vector_;
};

}//end of namespace schon

#endif /* MARGINALINFECTIONPROBABILITY_HPP_ */
