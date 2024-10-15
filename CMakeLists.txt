cmake_minimum_required(VERSION 3.25.1)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(MyProject LANGUAGES C CXX)

add_subdirectory(submodules/CMA-ES)
# install(DIRECTORY ${CMAKE_SOURCE_DIR}/submodules/CMA-ES/include/ DESTINATION include)
# # set_target_properties(cmaes PROPERTIES EXCLUDE_FROM_ALL TRUE)

# # Manually trigger the cmaes compilation when needed
# add_custom_target(build_cmaes
#     COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target cmaes
#     COMMENT "Building the cmaes library"
# )

# Options to be applied to all the projects
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

# Add subdirectories (add subprojects here, like your libraries and tests)
add_subdirectory(src/DPMathLib)

# Options to be applied to the dp_ libraries and the executables
    # Define the list of external libraries
    set(EXTERNAL_LIBS
        pthread 
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

add_subdirectory(src/dp_scalar2to2)
# add_subdirectory(src/dp_su2f2to2)
# add_subdirectory(src/dp_mssm2to2)

# We need to get the names of the libraries 
# To do so, we strip out everything before dp_ in the folder paths
set(MODELS_LIBS_NAMES "dp_scalar2to2")

# install(TARGETS DPMathLib ${MODELS_LIBS_NAMES} DESTINATION lib)

# Options to be applied only to the executables
# add_custom_target(install_libs
#     COMMAND cmake --build ${CMAKE_BINARY_DIR} --target install
#     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#     DEPENDS DPMathLib dp_scalar2to2
# )

add_subdirectory(src/Test)

# Options to be applied to the dp_ libraries and the executables
    # Define the list of external libraries
    set(EXTERNAL_LIBS
        pthread 
        ooptools 
        gfortran 
        gsl 
        gslcblas
        marty-lha
        cmaes
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


add_subdirectory(src/main)
# Assicurati che Test dipenda dall'installazione delle librerie
# add_dependencies(scalar_debug.x install_libs)

# Trova tutti i target definiti nel progetto Test
# get_property(TestTargets GLOBAL PROPERTY TARGETS)

# Filtro dei target che appartengono al progetto Test
# foreach(target IN LISTS TestTargets)
#     if("${target}" MATCHES "Test")  # Filtra i target con "Test" nel nome
#         add_dependencies(${target} install_libs)  # Aggiungi dipendenze
#     endif()
# endforeach()

# foreach(lib_path IN LISTS MODEL_LIB_NAMES)
#     # string(REGEX REPLACE ".*(dp_.*)$" "\\1" stripped_lib ${lib_path})
#     string(REGEX REPLACE ".*/(dp_.*)$" "src/\\1" stripped_lib ${lib_path})
#     list(APPEND MODELS_LIBS_NAMES ${stripped_lib})
# endforeach()