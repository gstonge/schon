#!/bin/bash
g++ -std=c++17 -o test_GroupSIS test_GroupSIS.cpp BipartiteNetwork.cpp GroupSIS.cpp Prevalence.cpp MarginalInfectionProbability.cpp -LSamplableSet/build/ -lsamplableset -ISamplableSet/ -g
