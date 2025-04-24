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

if [[ $(clang-format --version | cut -f 3 -d " ") == "19.1.7" ]]; then
    echo "Running format check"
    cmake --build ./build --target format-check || exit 1   
else
    echo "clang-format version is not 19.1.7, skipping format check"
fi

echo "Installing headers"
cmake --build ./build --target install_headers_all  || exit 1

echo "Building darkpack"
ninja -C build || exit 1
