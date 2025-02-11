#!/bin/bash

echo "Setting compilers"
export CC=gcc
export CXX=g++

echo "Testing read/write permissions"
mkdir build || exit 1

echo "Configuring cmake"
cmake -S . -B ./build \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=./build/.install \
      -DENABLE_CLANG_FORMAT=ON \
      -G Ninja || exit 1


echo "Running format check"
cmake --build ./build --target format-check || exit 1

echo "Installing headers"
cmake --build ./build --target install_headers_all  || exit 1

echo "Building documentation"
cmake --build ./build --target darkpack_doc || exit 1

echo "Building darkpack"
ninja -C build || exit 1
