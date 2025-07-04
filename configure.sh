#! /bin/bash

cmake -S . -B ./build \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=./build/.install \
      -DENABLE_CLANG_FORMAT=ON \
      -DENABLE_CLANG_TIDY=ON \
      -DENABLE_SUBMODULE_CMA=OFF \
      -G Ninja \
      -DMODELS_LIBS_NAMES="dp_u1f2to2"

if [[ $(clang-format --version | cut -f 3 -d " ") == "19.1.7" ]]; then
    cmake --build ./build --target format || exit 1   
else
    echo "clang-format version is not 19.1.7, skipping lining"
fi
# cmake --build ./build --target tidy || exit 1   
cmake --build ./build --target install_headers_all darkpack_doc 