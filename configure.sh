#! /bin/bash

libraries_list=(DPMathLib dp_scalar2to2)

cmake -S . -B ./build \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=./build/.install \
      -G Ninja

cd build

for elem in ${libraries_list[@]}
do
      cmake --build . --target install_headers_$elem
done