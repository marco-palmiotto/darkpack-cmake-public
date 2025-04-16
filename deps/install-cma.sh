#!/bin/sh

cd CMA-ES
prefix="$(pwd)/.install"
cmake -S . -B build
cd build
cmake --build . 
cmake --install . --prefix=${prefix}

