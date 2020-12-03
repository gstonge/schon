#!/bin/bash
g++ -std=c++17 -o test test_ComplexSIS.cpp BipartiteNetwork.cpp ComplexSIS.cpp -LSamplableSet/build/ -lsamplableset -ISamplableSet/
