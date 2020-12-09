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

#ifndef POWERLAWSIS
#define POWERLAWSIS

#include "ComplexSIS.hpp"

namespace schon
{//start of namespace schon

//class to simulate SIS process on networks
class PowerlawSIS : public ComplexSIS
{
public:
    //Constructor
    PowerlawSIS(const EdgeList& edge_list, double scale_recovery,
            double scale_infection, double shape_infection,
            const std::pair<double,double>& group_rate_bounds);
};

//constructor definiton
PowerlawSIS::PowerlawSIS(const EdgeList& edge_list, double scale_recovery,
            double scale_infection, double shape_infection,
            const std::pair<double,double>& group_rate_bounds) : ComplexSIS(
                edge_list, NULL, NULL, group_rate_bounds)
{
    //define the recovery/infection rates
    recovery_rate_ = [=](std::size_t n,std::size_t i) -> double
        {return scale_recovery*i;};
    infection_rate_ = [=](std::size_t n,std::size_t i) -> double
        {return scale_infection*(n-i)*pow(i,shape_infection);};
}

}//end of namespace schon

#endif /* POWERLAWSIS */

