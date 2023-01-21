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

#ifndef INFECTIOUS_SET_HPP_
#define INFECTIOUS_SET_HPP_

#include "Measure.hpp"
#include <vector>
#include <unordered_set>

namespace schon
{//start of namespace schon

class InfectiousSet : public Measure
{
public:
    //Constructor
    InfectiousSet();

    //Acessors
    std::vector<std::unordered_set<Node>> get_result() const;
    const std::string& get_name() const
        {return name_;}

    //Mutators
    void measure(ContagionProcess const * const ptr);
    void clear()
        {infectious_set_vector_.clear();}

private:
    //Members
    const std::string name_;
    std::vector<std::unordered_set<Node>> infectious_set_vector_;
};

}//end of namespace schon

#endif /* INFECTIOUS_SET_HPP_ */
