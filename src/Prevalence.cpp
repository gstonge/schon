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

#include "Prevalence.hpp"

using namespace std;

namespace schon
{//start of namespace schon

//constructor
Prevalence::Prevalence(size_t network_size): network_size_(network_size),
    count_(0), name_("prevalence")
{
}

//return the prevalence
double Prevalence::get_result() const
{
    return weight_/(count_*network_size_);
}

//perform a measure on the contagion process
void Prevalence::measure(
        ContagionProcess const * const ptr)
{
    weight_ += ptr->get_number_of_infected_nodes();
    count_ += 1;
}

}//end of namespace schon
