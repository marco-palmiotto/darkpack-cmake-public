#! /bin/bash

cmake -S . -B ./build \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=./build/.install \
      -DMODELS_LIBS_NAMES="dp_scalar2to2;dp_su2f2to2;dp_u1f2to2" \
      -DENABLE_CLANG_FORMAT=ON \
      -G Ninja

cmake --build ./build --target format && \
cmake --build ./build --target format-check && \
cmake --build ./build --target install_headers_all 
# cmake --build ./build --target darkpack_doc 