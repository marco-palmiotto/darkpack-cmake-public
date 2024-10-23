cmake_minimum_required(VERSION 3.10)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(MyProject LANGUAGES C CXX)

# Defining which models are to be compiled

macro(compile_all_libs)
# Call this macro to compile all the libraries
# We need to get the names of the libraries 
# To do so, we strip out everything before dp_ in the folder paths
    set(MODELS_LIBS_NAMES "")
    foreach(lib_path IN LISTS MODEL_LIB_NAMES)
        string(REGEX REPLACE ".*/(dp_.*)$" "src/\\1" stripped_lib ${lib_path})
        list(APPEND MODELS_LIBS_NAMES ${stripped_lib})
    endforeach()
endmacro()

set(MODELS_LIBS_NAMES "dp_scalar2to2")

# This submodule contains the CMA-ES algorithms for the parameter
# space sampling
add_subdirectory(submodules/CMA-ES)

# Options to be applied to ALL the projects
    # Set C++ standard globally
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    # Set C standard globally
    set(CMAKE_C_STANDARD 17)
    set(CMAKE_C_STANDARD_REQUIRED ON)

    # Set global compilation flags (applied to all targets)
    set(GENERIC_FLAGS "-Wall -Wextra -Wpedantic -Wno-deprecated-declarations -fPIC -fdiagnostics-color=always")
    set(MATH_OPTI_FLAGS "-fno-math-errno -ffinite-math-only -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fexcess-precision=fast -fno-signed-zeros -fno-trapping-math -fassociative-math -freciprocal-math")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GENERIC_FLAGS} ${MATH_OPTI_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GENERIC_FLAGS} ${MATH_OPTI_FLAGS} -Wshadow")

    # Debug mode flags (applied when building with Debug configuration)
    set(DEBUG_FLAGS "-g -O0") # Enable debug symbols, disable optimization
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")

    # Release mode flags (applied when building with Release configuration)
    set(RELEASE_FLAGS "-O3 -DNDEBUG") # Optimize for performance, disable debugging features
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${RELEASE_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")

    # Adding pthread to the external libs, to enable multithreading 
    set(EXTERNAL_LIBS pthread)

# Add subdirectories (add subprojects here, like your libraries and tests)
add_subdirectory(src/DPMathLib)

# Options to be applied to the dp_ libraries and the executables
    # Define the list of external libraries: 
    # Appending looptools, the fortran library, GSL, and MARTY
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} 
        ooptools 
        gfortran 
        gsl 
        gslcblas
        marty-lha
    )
    # Define the list of compile options (if needed)
    set(COMPILE_LINK_FLAGS
        -lpthread
        -looptools
        -lgfortran
        -lgsl
        -lgslcblas
        -lmarty-lha
    )

# Initiating the variable that contains all the libraries whose headers need to be installed before build
set(LISTS_COMMANDS_TO_INSTALL_HEADERS_ONLY DPMathLib)

foreach(dplib IN LISTS MODELS_LIBS_NAMES)
    # Adding the models to the projects
    add_subdirectory(src/${dplib}) 
    # Adding the command to install the header of the model to the dependencies of the global command defined next
    list(APPEND LISTS_COMMANDS_TO_INSTALL_HEADERS_ONLY install_headers_${dplib}) 
endforeach()

# Create a custom target for headers installation
add_custom_target(install_headers_all ALL DEPENDS ${LISTS_COMMANDS_TO_INSTALL_HEADERS_ONLY})

add_subdirectory(src/Test)

# Options to be applied to the dp_ libraries and the executables
# Appending CMA-ES to the list of external libs to be linked 
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} cmaes)
add_subdirectory(src/main)