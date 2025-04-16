INSTALLCMAPATH="$(pwd)/CMA-ES/build/.install"
INSTALLMARTYPATH="$(pwd)/marty-public/build/.install"

ALLPATHS=(${INSTALLCMAPATH} ${INSTALLMARTYPATH})

for single_path in ${ALLPATHS[@]}
do
  export PATH=$single_path/bin:$PATH
  export LD_LIBRARY_PATH=$single_path/lib
  export LIBRARY_PATH=$single_path/lib
  export C_INCLUDE_PATH=$single_path/include
  export CPATH=$single_path/include
done

