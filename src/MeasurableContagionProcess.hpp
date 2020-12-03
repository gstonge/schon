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

#ifndef MEASURABLECONTAGIONPROCESS_HPP_
#define MEASURABLECONTAGIONPROCESS_HPP_

#include <memory>
#include "ContagionProcess.hpp"
#include "Measure.hpp"
//#include "MarginalInfectionProbability.hpp"

namespace schon
{//start of namespace schon

//abstract class that allow measurement
class MeasurableContagionProcess : public ContagionProcess
{
public:
    //Constructor
    MeasurableContagionProcess(): measure_vector_() {};

    //Accessors
    const std::vector<std::shared_ptr<Measure>>& get_measure_vector() const
        {return measure_vector_;}

    //Mutators
    void add_measure(std::shared_ptr<Measure> ptr)
        {measure_vector_.push_back(ptr);}
    //void measure_marginal_infection_probability()
        //{add_measure(std::make_shared<MarginalInfectionProbability>(
                    //get_network().size()));}

protected:
    //Members
    std::vector<std::shared_ptr<Measure>> measure_vector_;
};

}//end of namespace schon

#endif /* MEASURABLECONTAGIONPROCESS_HPP_ */
