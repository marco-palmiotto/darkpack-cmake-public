#!/bin/bash

if [ $# -ne 2 ]
then 
  echo "This script requires two arguments:"
  echo " 1. The name of the numerical library"
  echo " 2. The name of the source file that calls MARTY, without the extension (.cpp)"
  exit 1
fi

echo Checking gcc version

# Get the version of GCC
gcc_version=$(gcc --version | grep -oP '\d+\.\d+\.\d+' | head -n1 | cut -f 1 -d '.')

echo You have installed $gcc_version

# Define the minimum required GCC version
required_version="12"

echo The required version is $required_version

# Compare the versions
if [ "$gcc_version" -ge "$required_version" ]; then
    # GCC version is greater than 12, export the variables
    export CC=gcc-12
    export CXX=g++-12
    export FC=gfortran-12

    echo Environmental variables have been updated
else
    export CXX=g++
    echo Your system is good
fi

# This script generates the numerical library named
libname=$1
# from the file 
source_mty=$2

[ ! -f $source_mty.cpp ] && echo "File $source_mty.cpp does not exist" && exit 2

# NOTE: to run this script, a configuration like the one in the commented lines 
#       has to be present in the .bashrc configuration file
#         # MARTY envirnmental variables
#         export INSTALLMARTYPATH=/home/$(USER)/installations/MARTY # path where MART is installed
#         export PATH=$INSTALLMARTYPATH/bin:$PATH # Find the command grafed
#         export LD_LIBRARY_PATH=$INSTALLMARTYPATH/lib:$LD_LIBRARY_PATH # Find libraries at run-time (Linux)
#         export LIBRARY_PATH=$INSTALLMARTYPATH/lib:$LIBRARY_PATH # Find libraries at compile-time
#         export C_INCLUDE_PATH=$INSTALLMARTYPATH/include:$C_INCLUDE_PATH # Find headers at compile time
#         export CPATH=$INSTALLMARTYPATH/include:$CPATH # Find headers at compile time

compile_MARTY_source()
{
  [ $# -ne 1 ] && echo "This function requires 1 argument : the name of the target file" && exit 1
  echo 'Generating C++ library via g++...'
  $CXX -std=c++17 -g -c $1.cpp -o $1.o -I$INSTALLMARTYPATH/include || exit 2
  [[ ! -f "$1.o" ]] && echo 'The file .o has not been generated. Compilation will stop now.' && exit 1
  echo 'Generating executable file linking with the libraries of MARTY...'
  $CXX -std=c++17 -g -o $1.x $1.o -L$INSTALLMARTYPATH/lib -lmarty || exit 2
  [[ ! -f "$1.x" ]] && echo 'The file .x has not been generated. Compilation will stop now.' && exit 1
}

[ -e "$libname"_graphs ] || mkdir "$libname"_graphs
[ -e "$libname"_graphs/README ] || $(echo "This directory contains the Feynman diagrams of the processes computed in $source_mty.cpp" > "$libname"_graphs/README)

# Deleting the file with the calcHep-SuperIso correspondance
# It will be opened in append mode next
rm -f "data/processes_chep.txt"

# We need those folders to create files in them while $source_mty.x is run
[ -d $libname/include ] || mkdir -p $libname/include
[ -d $libname/src ] || mkdir -p $libname/src

rm "$source_mty.x" "$source_mty.o"
compile_MARTY_source $source_mty

./$source_mty.x || exit 3


echo "Library configuration"

auxfolder=auxiliary_library
auxfolder_model=$auxfolder/$libname
exit_code=0

print_how()
{
  echo "This script reads needs 1 argument(s):"
  echo "   1. the first is name of the library"
}

generate_param_t_comparison()
{
# This function has to be launched inside the $auxfolder_model
  local libpath=../../../$libname/include/params.h

  grep "csl::InitSanitizer<[[:alnum:]_]*> [[:alnum:]_]* { \"[[:alnum:]_]*\" };" $libpath | cut -f 2 -d '>' | cut -f 1 -d '{' | sed "s/^ \(.*\)$/if(!(p1.\1==p2.\1))return false;/g" > central_test.txt
  
  cat ../blank_files/compare_param_t.head central_test.txt > compare.cpp

  echo 'return true;}}' >> compare.cpp
  
  rm central_test.txt

}

generate_copy_to_c_struct()
{
# This function has to be launched inside the $auxfolder_model
  local inputfile=../../../$libname/include/params.h
  local outputfile=copy_to_c_struct.hpp
  
  echo Determining sij

  local ilist=( )
  local jlist=( )

  for i in $(seq 1 4)
  do
    for j in $(seq 1 4)
    do
      local filtered=$(grep "csl::InitSanitizer<real_t> s_$i$j" $inputfile )
      if [[ $filtered != '' ]]
      then 
        ilist+=($i)
        jlist+=($j)
        echo Found s_$i$j
#         echo Updated lists are, for the is
#         for elem in ${ilist[@]} 
#         do 
#           echo --$elem--
#         done
#         echo For the js
#         for elem in ${jlist[@]}
#         do 
#           echo --$elem--
#         done
      fi
    done
  done

  local npairs=${#ilist[@]}
  echo "Found $npairs pairs" 

#   echo Writing the first function 

  echo "namespace $libname" > $outputfile
  echo "{" >> $outputfile
  echo "inline void update_kinematics(Param_t &input, const double sij[5][5])" >> $outputfile
  echo "{" >> $outputfile

  npairs=$(($npairs - 1))

  for counter in $(seq 0 $npairs)
  do
    echo "  input.s_${ilist[$counter]}${jlist[$counter]} = sij[${ilist[$counter]}][${jlist[$counter]}];" >> $outputfile
  done

  echo "}" >> $outputfile

  echo "inline cparam_s copy_to_c_struct(Param_t &input)">> $outputfile
  echo "{">> $outputfile
  echo "  cparam_s output;">> $outputfile
  echo "#define CONVERT_COMPLEX(X) ( (static_cast<std::complex<double>>(X)).real() + _mty_I*(static_cast<std::complex<double>>(X)).imag() )" >> $outputfile
  ##########

  local realparameters=()
  local complexparameters=()

  realparameters=$(grep "  csl::InitSanitizer<real_t>" $inputfile | cut -f 2 -d '>' | cut -f 1 -d '{' | sed 's/ //g' | tr '\n' ' ')

  complexparameters=$(grep "  csl::InitSanitizer<complex_t>" $inputfile | cut -f 2 -d '>' | cut -f 1 -d '{' | sed 's/ //g' | tr '\n' ' ')

  for param in ${realparameters[@]}
  do 
    echo "  output.$param = input.$param;" >> $outputfile
  done

  for param in ${complexparameters[@]}
  do 
    echo "  output.$param = CONVERT_COMPLEX(input.$param);" >> $outputfile
  done

  #########

  echo "#undef CONVERT_COMPLEX">> $outputfile
  echo "  return output;">> $outputfile
  echo "}">> $outputfile

  echo "inline cparam_s update_c_kinematics(const Param_t &input, cparam_s output)">> $outputfile
  echo "{">> $outputfile

  for counter in $(seq 0 $npairs)
  do
    echo "  output.s_${ilist[$counter]}${jlist[$counter]} = input.s_${ilist[$counter]}${jlist[$counter]};" >> $outputfile
  done

  echo "  return output;" >> $outputfile

  echo "}" >> $outputfile

  echo "} // End of namespace $libname" >> $outputfile
}



[ -d "$auxfolder_model" ] || mkdir -p "$auxfolder_model"

echo "Generating compare.cpp and copy_to_c_struct.hpp"
curdir=$(pwd)
cd $auxfolder_model
generate_param_t_comparison 
generate_copy_to_c_struct
cd $curdir
cp $auxfolder/blank_files/compare.hpp $auxfolder_model/compare.hpp

