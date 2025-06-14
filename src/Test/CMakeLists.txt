cmake_minimum_required(VERSION 3.10)

function(print)
    foreach(var ${ARGN})
        message("${var} = ${${var}}")
    endforeach(var ${ARGN})
endfunction()

project(test)

file(GLOB MODEL_LIB_NAMES ${PROJECT_SOURCE_DIR}/../dp_* )
set(INVARIANT_LIB_NAME DPMathLib)

file(GLOB SOURCE_FILES_LIST ${PROJECT_SOURCE_DIR}/*.cpp)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    # Independent compilation

    # Define the list of external libraries
    set(EXTERNAL_LIBS
                        pthread 
                        ooptools 
                        gfortran 
                        gsl 
                        gslcblas
        )
    # Define the list of compile options (if needed)
    set(COMPILE_LINK_FLAGS
                        -lpthread
                        -looptools
                        -lgfortran
                        -lgsl
                        -lgslcblas
        )
    
    # Flags definition

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

    # We need to get the names of the libraries to link to 
    # the executables.
    # To do so, we strip out everything before dp_ in the folder paths
    set(MODELS_LIBS_NAMES "")
    foreach(lib_path IN LISTS MODEL_LIB_NAMES)
        string(REGEX REPLACE ".*(dp_.*)$" "\\1" stripped_lib ${lib_path})
        list(APPEND MODELS_LIBS_NAMES ${stripped_lib})
    endforeach()
endif()

print(MODEL_LIB_NAMES MODELS_LIBS_NAMES SOURCE_FILES_LIST)

# 1. Generate an executable for each .cpp source file with a .x extension
foreach(source_file IN LISTS SOURCE_FILES_LIST)
    # Get the filename without the path and the extension
    get_filename_component(executable_name ${source_file} NAME_WE)
    # Add the executable with .x extension
    add_executable(${executable_name}.x ${source_file})

    message("${executable_name}.x ${source_file}")

    # Include directories and link libraries for the executable
    target_include_directories(${executable_name}.x PUBLIC
        ${CMAKE_INSTALL_PREFIX}/include
    )
    target_link_libraries(${executable_name}.x ${INVARIANT_LIB_NAME} ${EXTERNAL_LIBS} ${MODELS_LIBS_NAMES})
endforeach()