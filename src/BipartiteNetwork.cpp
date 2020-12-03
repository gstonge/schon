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

#include "BipartiteNetwork.hpp"
#include <numeric>

using namespace std;

namespace schon
{//start of namespace schon

//Constructor of the class provided an edge list
BipartiteNetwork::BipartiteNetwork(const EdgeList& edge_list) :
	node_adjacency_list_(), group_adjacency_list_(), nodes_(), groups_(),
    min_membership_(0), max_membership_(0), min_group_size_(0),
    max_group_size_(0)
{
	size_t nb_nodes = 0;
	size_t nb_groups = 0;
	//Determine the number of nodes and groups (brut force)
	for (size_t i = 0; i < edge_list.size(); i++)
    {
    	if (edge_list[i].first > nb_nodes)
    	{
    		nb_nodes = edge_list[i].first;
    	}
    	if (edge_list[i].second > nb_groups)
    	{
    		nb_groups = edge_list[i].second;
    	}
    }
    nb_nodes += 1; //the label starts to 0 by convention
    nb_groups += 1; //the label starts to 0 by convention

    //Initialize adjacency lists, nodes an groups
    node_adjacency_list_ = NodeAdjacencyList(nb_nodes, vector<Group>());
    group_adjacency_list_ = GroupAdjacencyList(nb_groups, vector<Node>());
    nodes_ = vector<Node>(nb_nodes);
    groups_ = vector<Node>(nb_groups);
    iota(nodes_.begin(),nodes_.end(),0);
    iota(groups_.begin(),groups_.end(),0);

    for (auto & edge : edge_list)
    {
        node_adjacency_list_[edge.first].push_back(edge.second);
        group_adjacency_list_[edge.second].push_back(edge.first);
    }

    //Determine min and max membership
    for (Node node : nodes_)
    {
        if (node == 0)
        {
            min_membership_ = membership(node);
        }
        if (membership(node) < min_membership_)
        {
            min_membership_ = membership(node);
        }
        if (membership(node) > max_membership_)
        {
            max_membership_ = membership(node);
        }
    }

    //Determine min and max  group size
    for (Group group : groups_)
    {
        if (group == 0)
        {
            min_group_size_ = group_size(group);
        }
        if (group_size(group) < min_group_size_)
        {
            min_group_size_ = group_size(group);
        }
        if (group_size(group) > max_group_size_)
        {
            max_group_size_ = group_size(group);
        }
    }
}

}//end of namespace schon
