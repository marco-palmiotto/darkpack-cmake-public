cmake_minimum_required(VERSION 3.10)

# setting the namespace for the specialised library
set(__SPEC_LIB_NAME__ scalar2to2)

function(print)
    foreach(var ${ARGN})
        message("${var} = ${${var}}")
    endforeach(var ${ARGN})
endfunction()

project(dp_${__SPEC_LIB_NAME__} VERSION 1.0 LANGUAGES C CXX)
set(EXEC_MAIN_NAME main_final.x)

set(INVARIANT_LIB_NAME DPMathLib)

set(SPEC_LIB_PATH ${PROJECT_SOURCE_DIR}/../${__SPEC_LIB_NAME__})
set(TEMPLATED_LIB_PATH ${PROJECT_SOURCE_DIR}/../TemplatedLib)

# The invariant is not compiled here, but as standalone
# It needs to be linked for the executables

set(INVARIANT_LIB_HEADER_PATH ${PROJECT_SOURCE_DIR}/../${INVARIANT_LIB_NAME})    
set(INVARIANT_LIB_PATH ${INVARIANT_LIB_HEADER_PATH}/build/lib) 

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
# Set global compilation flags (applied to all targets)
set(GENERIC_FLAGS "-Wall -Wextra -Wpedantic -Wno-deprecated-declarations -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GENERIC_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GENERIC_FLAGS}")
# Debug mode flags (applied when building with Debug configuration)
set(DEBUG_FLAGS "-g -O0") # Enable debug symbols, disable optimization
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
# Release mode flags (applied when building with Release configuration)
set(RELEASE_FLAGS "-O3 -DNDEBUG") # Optimize for performance, disable debugging features
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${RELEASE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")

# Creating configuration file with library-specific macros
# Note that the destination header file must have a generic name,
# in order for us to write the "templated library" in a way
# independent of the "specialization library"
configure_file(${PROJECT_SOURCE_DIR}/config.hpp.in 
               ${PROJECT_SOURCE_DIR}/config.hpp
               )

set(SOURCE_LIB_DIR ${PROJECT_SOURCE_DIR})

print(CMAKE_BINARY_DIR SPEC_LIB_PATH TEMPLATED_LIB_PATH INVARIANT_LIB_HEADER_PATH SOURCE_LIB_DIR)

# Automatically gather all .cpp files in all projects to be compiled at once
file(GLOB LIB_CXX_SOURCES ${SPEC_LIB_PATH}/src/*.cpp 
                           ${TEMPLATED_LIB_PATH}/*.cpp 
                           ${SOURCE_LIB_DIR}/*.cpp)

# Automatically gather all .c files in all projects to be compiled at once
file(GLOB LIB_C_SOURCES ${SPEC_LIB_PATH}/clib/*.c
                        ${TEMPLATED_LIB_PATH}/clib/*.c
                        )

# Automatically gather all the CXX header files in all projects to be installed later on
file(GLOB LIB_CXX_HEADERS  ${SPEC_LIB_PATH}/include/*.h
                           ${SPEC_LIB_PATH}/include/*.hpp
                           ${TEMPLATED_LIB_PATH}/*.h
                           ${TEMPLATED_LIB_PATH}/*.hpp
                           ${CMAKE_CURRENT_SOURCE_DIR}/*.h
                           ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)

# Automatically gather all the C header files in all projects to be installed later on
file(GLOB LIB_C_HEADERS ${SPEC_LIB_PATH}/clib/*.h
                        ${TEMPLATED_LIB_PATH}/clib/*.h)

# print(LIB_CXX_SOURCES LIB_C_SOURCES LIB_CXX_HEADERS LIB_C_HEADERS)      

# For debugging: gathering all the header files of the invariant library
file(GLOB INVARIANT_LIB_HEADER_FILES ${INVARIANT_LIB_HEADER_PATH}/*.h
                                     ${INVARIANT_LIB_HEADER_PATH}/*.hpp)

# print(INVARIANT_LIB_HEADER_FILES)
            
# Defining all the targets:
    # Add the library and include all the .cpp files
add_library(${PROJECT_NAME} ${LIB_CXX_SOURCES} ${LIB_C_SOURCES})
    # Add the executable
add_executable(${EXEC_MAIN_NAME} script/main.cpp)

# List of targets in pure C++
set(TARGETS ${PROJECT_NAME} ${EXEC_MAIN_NAME})

# For each CXX target: 
#   - specify the include directoryu
#   - link with the invariant library
foreach(target IN LISTS TARGETS)
    target_include_directories(${target}  PUBLIC
        ${INVARIANT_LIB_HEADER_PATH}         # C/C++ lib files of the invariant library 
        ${SPEC_LIB_PATH}/include             # C++ lib files of the specialised library 
        ${SPEC_LIB_PATH}/clib                # C   lib files of the specialised library 
        ${TEMPLATED_LIB_PATH}                # C/C++ lib files of the templated library
        ${CMAKE_CURRENT_SOURCE_DIR}          # C/C++ lib files
    )
    target_link_directories(${target}  PUBLIC
        ${INVARIANT_LIB_PATH}
    )
    target_compile_options(${target} PRIVATE 
        ${COMPILE_LINK_FLAGS}
    )
endforeach()

# for the library-targets
#    Specify output folder 
set_target_properties(${PROJECT_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)
#    Linking external library
target_link_libraries(${PROJECT_NAME} PRIVATE
    ${EXTERNAL_LIBS}
    ${INVARIANT_LIB_NAME}
)

# Copying the headers of the C++ library
install(FILES ${LIB_CXX_HEADERS} DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
# Copying the headers of the C library
install(FILES ${LIB_C_HEADERS} DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}/C)

# For the executables
#    Manual link to the path of the libraries
target_link_directories(${EXEC_MAIN_NAME}  PRIVATE ${CMAKE_BINARY_DIR}/lib/)
#    Manual link of the library
target_link_libraries(${EXEC_MAIN_NAME} PRIVATE
    ${PROJECT_NAME} 
    ${EXTERNAL_LIBS}
    ${INVARIANT_LIB_NAME}
)


