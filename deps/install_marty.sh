cd marty-public 
mkdir build
cd build 
cmake .. -DCMAKE_INSTALL_PREFIX=.install -G Ninja 
cmake --build . 
cmake --install .