#! /usr/bin/bash
g++ verlet.cpp schrodinger.cpp bisection.cpp
./a.out > $1
