#! /bin/bash

cmake -S . -B ./build \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=./build/.install \
      -DENABLE_CLANG_FORMAT=ON \
      -DENABLE_SUBMODULE_CMA=OFF \
      -G Ninja \
      -DMODELS_LIBS_NAMES="dp_u1f2to2"

cmake --build ./build --target format && \
cmake --build ./build --target format-check install_headers_all darkpack_doc #&& \
# cmake --build ./build --target install_headers_all 
# cmake --build ./build --target darkpack_doc 