# Changelog

## 21st february 2024 (v1.2)

Main changes:

- Upgraded `DarkPACK` to support the latest version of `MARTY`

- Overloaded `<<` for `Param_t`

- `MARTY` itself has been upgraded to support explicitly the assignment operator of `InitSanitizer`, by ensuring that its name remains untouched, allowing for consistent debugging when an exception is thrown

# `CMakeLists.txt` Documentation

## Overview

This `CMakeLists.txt` file is used to configure and build the `MyProject` project. It sets up the necessary build options, finds required packages, and defines custom targets for documentation generation and code formatting.

## Functions

### `print`

Prints the values of the given variables.

### `list_all_models`

Lists all models in the `src` directory that start with `dp_`.

### `should_exclude`

Checks if a file should be excluded based on the defined patterns.

## Options

### `BUILD_DOCUMENTATION`

Option to enable or disable documentation generation with Doxygen.

### `ENABLE_SUBMODULE_CMA`

Option to enable or disable the compilation of the CMA submodule and dependent libraries.

### `ENABLE_CLANG_FORMAT`

Option to enable or disable clang-format to modify files.

## Variables

### `MODELS_LIBS_NAMES`

List of model libraries to compile.

### `EXTERNAL_LIBS`

List of external libraries to be linked.

### `COMPILE_LINK_FLAGS`

List of compile options for linking external libraries.

### `LISTS_COMMANDS_TO_INSTALL_HEADERS_ONLY`

List of commands to install headers only.

### `EXCLUDE_PATTERNS`

Patterns to exclude from clang-format.

## Targets

### `install_headers_all`

Custom target for headers installation.

### `darkpack_doc`

Custom target for documentation generation with Doxygen.

### `format`

Custom target to run clang-format to change files.

### `format-check`

Custom target to run clang-format to check code style compliance.

## Subdirectories

### `src/DPMathLib`

Adds the `DPMathLib` subdirectory.

### `src/Test`

Adds the `Test` subdirectory.

### `submodules/CMA-ES`

Adds the `CMA-ES` submodule if `ENABLE_SUBMODULE_CMA` is enabled.

## Documentation Generation

If `BUILD_DOCUMENTATION` is enabled and Doxygen is found, a custom target for documentation generation is added.

## Clang-Format

If `ENABLE_CLANG_FORMAT` is enabled, custom targets for running clang-format are added.

## Conclusion

This `CMakeLists.txt` file sets up the build environment for `MyProject`, including options for documentation generation, code formatting, and submodule compilation. It defines custom functions, options, variables, and targets to streamline the build process.


# Requirements 

In order to use this package, `MARTY` has to be installed and configured on your machine. 
No further external dependencies are required.

We encourage users to read carefully this section on requirements before proceeding with the installation. It may save you some time, especially if you are on Linux and you are installing `MARTY` for the first time.


You can copy this package in the directory that you prefer:

	wget https://gitlab.in2p3.fr/darkpack/darkpack-public/-/archive/master/darkpack-public-master.tar.gz
	tar -xf darkpack-public-master.tar.gz
	cd darkpack-public-master
On Linux and Mac it is possible to download, install and configure `MARTY` by running the (unofficial) script `install_marty.sh` in the folder where you wish to install `MARTY`. 
The script will extract and compile the version on `MARTY`
contained in the file `marty-public.tar.gz`,
provided that the dependencies are installed.
You can find all the details about installations and dependencies [on MARTY's website](https://marty.in2p3.fr).

If you want to install the latest version, you can open the script and modify it. There is already a function to do it. For any problem with the installation of `MARTY` we recommend to visit [MARTY's website](https://marty.in2p3.fr). 

The following procedure refers to the `mssm2to2` library. The same results can be obtained for any library (named e.g. `newlib`), simply replacing `mssm2to2` with `newlib` and vice versa in what follows.

Note that all the scripts are written for a `bash` environment. If you wish to use another shell, you may need to adjust the scripts. It is important that the environmental variable `INSTALLMARTYPATH` is defined and that is the path where `MARTY` has been built. 
If you used the script to compile `MARTY`, it will help you to set this variable.

Once `MARTY` is installed, you are ready to use this package!

The manual can be found in this repository and also [here](https://arxiv.org/abs/2211.10376).

# Setup of the library

## Generation of new models
In order to generate a new model, you need to 
1. place its source file in the folder `src/models`;
2. compile it, by calling 

        ./lib_generate.sh <1> <2>
    
    where `<1>` is the name of the library you want as output, and `<2>` is the name of the `.cpp` file (without its extension)
3. now the model is generated, you need to copy it into `src`. To do so, run the script

        ./lib_copy.sh <1>
    where `<1>` is the same as before. This script will also genrate the `CMakeLists.txt` needed in the projects.

## Compiling with the new-generated model

4. you need to add your new model into the `CMakeLists.txt` in the root folder. Add the line:

        add_subdirectory(src/dp_<1>)
     where `<1>` is the same as before, after the other projects are added.
5. you can compile all the models together by using `cmake`. Otherwise, you can use the scripts we provide

        ./configure.sh
        ./build.sh
    which will, respectively, configure and build `DarkPACK` in release mode, and using `ninja` as build system.