#! /bin/bash

cmake -S . -B ./build \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=./build/.install \
      -G Ninja
      
cmake --build ./build --target install_headers_all
