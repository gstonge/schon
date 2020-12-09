#!/bin/bash
g++ -std=c++17 -o test_ComplexSIS test_ComplexSIS.cpp BipartiteNetwork.cpp ComplexSIS.cpp Prevalence.cpp MarginalInfectionProbability.cpp -LSamplableSet/build/ -lsamplableset -ISamplableSet/ -g
