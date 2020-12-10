import matplotlib.pyplot as plt
from _schon import PowerlawGroupSIS

#structure : all individuals belong to two groups
N = 500
edge_list = []
for node in range(N):
    edge_list.append((node,0))
    edge_list.append((node,1))

#infection parameter
recovery_rate = 1.
scale_infection = 1.2*10**(-3)
shape_infection = 1. #linear
rate_bounds = (10**(-3),100)
initial_infected_fraction = 0.05
seed = 42
nb_history = 50

cont = PowerlawGroupSIS(edge_list,recovery_rate,scale_infection,
                        shape_infection,rate_bounds)
cont.infect_fraction(initial_infected_fraction)
# cont.seed(seed) #optional
cont.initialize_history(nb_history)

#define some measures
cont.measure_prevalence()
cont.measure_marginal_infection_probability()

#evolve in the quasistationary state without measuring (burn-in)
dt = 100
dec_dt = 0.5
cont.evolve(dt,dec_dt,measure=False,quasistationary=True)

#evolve and measure
dt = 1000
dec_dt = 1
cont.evolve(dt,dec_dt,measure=True,quasistationary=True)

#print the result measure
for measure in cont.get_measure_vector():
    name = measure.get_name()
    if name == "prevalence":
        print("----------------")
        print(name)
        print("----------------")
        print(measure.get_result())
    elif name == "marginal_infection_probability":
        plt.hist(measure.get_result())
        plt.show()


