#!/bin/bash

[ -d build/submodules ] &&  mv build/submodules temp
rm -rf build
mkdir build 
mv temp build/submodules