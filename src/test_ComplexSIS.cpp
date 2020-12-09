#include "PowerlawSIS.hpp"
#include <iostream>

using namespace std;
using namespace schon;

int main(int argc, const char *argv[])
{
    //make one big group
    int n = 1000;
    double scale_rec = 1.;
    double scale_inf = 0.001;
    double shape_inf = 1.;
    double fraction = 0.5;
    EdgeList edge_list;
    for( int j = 0; j < n; j++)
    {
        edge_list.push_back(make_pair(j,0));
    }
    PowerlawSIS cont(edge_list,scale_rec,scale_inf,shape_inf,
            make_pair(1.,scale_inf*(n*n/4)+n*scale_rec));

    //measure prevalence and marginal infection prob
    cont.measure_prevalence();
    cont.measure_marginal_infection_probability();

    //infect a fraction of the nodes
    cont.infect_fraction(fraction);
    cont.initialize_history();

    //make it evolve for a time, then evolve and measure
    double dt1 = 200.;
    double dt2 = 20.;
    cont.evolve(dt1,0.5,false,true);
    cont.evolve(dt2,1.,true,true);
    vector<shared_ptr<Measure>> measure_vector = cont.get_measure_vector();
    shared_ptr<Prevalence> prevalence_ptr = dynamic_pointer_cast<Prevalence>(
            measure_vector[0]);
    double prevalence = prevalence_ptr-> get_result();
    cout << prevalence << endl;


    return 0;
}
