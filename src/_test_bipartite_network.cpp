#include "BipartiteNetwork.hpp"
#include <iostream>

using namespace std;
using namespace schon;

int main(int argc, const char *argv[])
{
    EdgeList edge_list;
    edge_list.push_back(make_pair(0,0));
    edge_list.push_back(make_pair(0,1));
    edge_list.push_back(make_pair(1,1));
    edge_list.push_back(make_pair(0,2));
    BipartiteNetwork net(edge_list);
    cout << net.membership(0) << endl;
    cout << net.group_size(1) << endl;
    cout << net.min_membership() << endl;
    cout << net.max_membership() << endl;
    cout << net.min_group_size() << endl;
    cout << net.max_group_size() << endl;


    return 0;
}
