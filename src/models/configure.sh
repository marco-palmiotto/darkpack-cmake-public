#! /bin/bash

# Library configuration

auxfolder=auxiliary_library

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
  
  cat compare_param_t.head central_test.txt > compare.cpp

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

  echo "} // End of namespace __SPEC_LIB_NAME__" >> $outputfile
}

# Assigning the name for the numerical library
if [ $# -ge 1 -a -d ../$1 ]
then 
  libname=$1 
else
  echo "You need to specify a valid name for the library"
  print_how
  exit 2
fi

auxfolder_model=$auxfolder/$libname
exit_code=0

[ -d "$auxfolder_model" ] || mkdir -p "$auxfolder_model"

echo "Generating compare.cpp and copy_to_c_struct.h"
curdir=$(pwd)
cd $auxfolder_model
generate_param_t_comparison 
generate_copy_to_c_struct
cd $curdir
