#include "ComplexSIS.hpp"
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

    return 0;
}
