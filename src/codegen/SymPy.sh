#! /bin/bash
g++ SympyTest.cpp -o Sympy -lm -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/unsupported/ -I ~/lib/
./Sympy
