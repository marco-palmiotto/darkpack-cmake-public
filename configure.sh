#! /bin/bash

cmake -S . -B ./build \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=./build/.install \
      -G Ninja

cd build

cmake --build . --target install_headers_all
