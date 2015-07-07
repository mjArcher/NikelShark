#!/bin/bash 
g++ LibraryTest.cpp SquareTensor3.cpp ElasticPrimState.cpp -I ~/Libraries/ -I ~/Libraries/eigen/ -I ~/Libraries/eigen-dev/ -I ~/Libraries/eigen-dev/unsupported/ -lblitz -o test
./test

