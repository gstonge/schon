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

#ifndef BIPARTITENETWORK_HPP_
#define BIPARTITENETWORK_HPP_

#include <utility>
#include <vector>
#include <memory>

namespace schon
{//start of namespace schon

typedef unsigned int Node;
typedef unsigned int Group;
typedef std::vector<std::pair<Node,Group>> EdgeList;
typedef std::vector<std::vector<Group> > NodeAdjacencyList;
typedef std::vector<std::vector<Node> > GroupAdjacencyList;


//Structure representing an undirected network
class BipartiteNetwork
{
public:
    //Constructor
    BipartiteNetwork(const EdgeList& edge_list);

    //Accessors
    std::size_t min_membership() const
        {return min_membership_;}
    std::size_t max_membership() const
        {return max_membership_;}
    std::size_t min_group_size() const
        {return min_group_size_;}
    std::size_t max_group_size() const
        {return max_group_size_;}

    std::size_t membership(Node node) const
    	{return node_adjacency_list_[node].size();}
    std::size_t group_size(Group group) const
    	{return group_adjacency_list_[group].size();}

    std::size_t size() const
        {return node_adjacency_list_.size();}
    std::size_t number_of_nodes() const
        {return node_adjacency_list_.size();}
    std::size_t number_of_groups() const
        {return group_adjacency_list_.size();}

    const std::vector<Node>& group_members(Group group) const
    	{return group_adjacency_list_[group];}
    const std::vector<Group>& adjacent_groups(Node node) const
    	{return node_adjacency_list_[node];}
    const std::vector<Node>& nodes() const
        {return nodes_;}
    const std::vector<Group>& groups() const
        {return groups_;}

private:
    //Members
    NodeAdjacencyList node_adjacency_list_;
    GroupAdjacencyList group_adjacency_list_;
    std::vector<Node> nodes_;
    std::vector<Group> groups_;
    std::size_t min_membership_;
    std::size_t max_membership_;
    std::size_t min_group_size_;
    std::size_t max_group_size_;

};

}//end of namespace schon

#endif /* BIPARTITENETWORK_HPP_ */
