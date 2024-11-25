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

# file list of the files to copy unconditionally
file_list=("config.hpp.in" ) 
for filename in ${file_list[@]}
do
    echo Copying $filename in $darkpack_lib_destination_path
    cp auxiliary_library/blank_files/$filename $darkpack_lib_destination_path/$filename
done

# file to be copied only if not present already
filename=leshouchesfrommarty.hpp
echo Checking if $filename is already present in the library
if [ -f $darkpack_lib_destination_path/$filename ]
then
    echo $darkpack_lib_destination_path/$filename already exists
else
    echo $darkpack_lib_destination_path/$filename does not exist, copying the template file
    cp auxiliary_library/blank_files/$filename $darkpack_lib_destination_path/$filename
fi
