#!/bin/bash

if [ $# -ne 1 ]
then 
  echo "This script requires one argument:"
  echo " 1. The name of the numerical library"
  exit 1
fi

[ ! -d $1 ] && echo "The library does not exist" && exit 2

barebone_lib_destination_path=../$1
if [ ! -d $barebone_lib_destination_path ]
then
    echo Creating the new library
    cp -r $1 $barebone_lib_destination_path
else
    echo A library already exists in /src, skipping.
    echo NOTE: you need to run 
    echo rm -rf $barebone_lib_destination_path
    echo if you want to delete it and replace it with the new one
fi

echo Generating CMakeLists.txt in $barebone_lib_destination_path
sed "s/bsm2to2/$1/g" auxiliary_library/blank_files/CMakeLists.txt > $barebone_lib_destination_path/CMakeLists.txt

darkpack_lib_destination_path=../dp_$1
if [ ! -d $darkpack_lib_destination_path ]
then 
    echo Creting the new DarkPACK library
    mkdir $darkpack_lib_destination_path
fi

dp_files=$(echo auxiliary_library/$1/*)
for filename in ${dp_files[@]}
do
    echo Processing $filename
    cp $filename $darkpack_lib_destination_path/
done

echo Generating CMakeLists.txt in $darkpack_lib_destination_path
sed "s/bsm2to2/$1/g" auxiliary_library/blank_files/dp/CMakeLists.txt > $darkpack_lib_destination_path/CMakeLists.txt