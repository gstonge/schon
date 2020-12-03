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

#ifndef CONTAGIONPROCESS_HPP_
#define CONTAGIONPROCESS_HPP_

#include <vector>
#include "BipartiteNetwork.hpp"
#include <unordered_set>
#include <unordered_map>
#include "utility.hpp"

namespace schon
{//start of namespace schon

enum NodeState {S, I, Count};
const unsigned int STATECOUNT = static_cast<unsigned int>(NodeState::Count)

typedef std::vector<std::vector<Node>> GroupState; //NodeState is entry
typedef std::unordered_map<Node,std::size_t> GroupStatePosition;

//abstract class with minimal structure for the simulation of contagions
class ContagionProcess
{
public:
    //Constructor
    ContagionProcess() {};

    //Accessors
    virtual const std::vector<NodeState>& get_node_state_vector() const = 0;
    virtual std::size_t get_number_of_infected_nodes() const = 0;
    virtual const std::unordered_set<Node>& get_infected_node_set(
            ) const = 0;
    virtual double get_lifetime() const = 0;
};

}//end of namespace schon

#endif /* CONTAGIONPROCESS_HPP_ */
