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

#include "GroupSIS.hpp"

namespace schon
{//start of namespace schon

//class to simulate SIS process on networks
class PowerlawGroupSIS : public GroupSIS
{
public:
    //Constructor
    PowerlawGroupSIS(const EdgeList& edge_list, double recovery_rate,
            double scale_infection, double shape_infection,
            const std::pair<double,double>& rate_bounds);
};

//constructor definiton
PowerlawGroupSIS::PowerlawGroupSIS(const EdgeList& edge_list, double recovery_rate,
            double scale_infection, double shape_infection,
            const std::pair<double,double>& rate_bounds) : GroupSIS(
                edge_list, recovery_rate, NULL, rate_bounds)
{
    //define the recovery/infection rates
    infection_rate_ = [=](std::size_t n,std::size_t i) -> double
        {return scale_infection*(n-i)*pow(i,shape_infection);};
}

}//end of namespace schon

#endif /* POWERLAWSIS */

