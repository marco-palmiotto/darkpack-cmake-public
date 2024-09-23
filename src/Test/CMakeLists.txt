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

# Strip out everything before dp_ in the folder paths
set(STRIPPED_MODEL_LIB_NAMES "")
foreach(lib_path IN LISTS MODEL_LIB_NAMES)
    string(REGEX REPLACE ".*(dp_.*)$" "\\1" stripped_lib ${lib_path})
    list(APPEND STRIPPED_MODEL_LIB_NAMES ${stripped_lib})
endforeach()

print(MODEL_LIB_NAMES STRIPPED_MODEL_LIB_NAMES SOURCE_FILES_LIST)

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

# Using the installed path for the invariant library's header
set(INVARIANT_LIB_HEADER_PATH ${CMAKE_BINARY_DIR}/include/${INVARIANT_LIB_NAME})  

# 1. Generate an executable for each .cpp source file with a .x extension
foreach(source_file IN LISTS SOURCE_FILES_LIST)
    # Get the filename without the path and the extension
    get_filename_component(executable_name ${source_file} NAME_WE)
    # Add the executable with .x extension
    add_executable(${executable_name}.x ${source_file})

    message("${executable_name}.x ${source_file}")

    # Include directories and link libraries for the executable
    target_include_directories(${executable_name}.x PUBLIC
        ${INVARIANT_LIB_HEADER_PATH}
        ${CMAKE_BINARY_DIR}/include
    )
    target_link_libraries(${executable_name}.x ${INVARIANT_LIB_NAME} ${EXTERNAL_LIBS} ${STRIPPED_MODEL_LIB_NAMES})
endforeach()

# add_executable(${EXEC_NAME_0} main.cpp)

# # for each executable, we perform some operations:
# # 1. defining new header location
# target_include_directories(${EXEC_NAME_0} PUBLIC
#     ${INVARIANT_LIB_HEADER_PATH}
#     ${CMAKE_BINARY_DIR}/include
# )

# # 2. manual link to the path of the libraries
# target_link_directories(${EXEC_NAME_0} PRIVATE ${CMAKE_BINARY_DIR}/lib/)
# # 3. manual link of the library
# target_link_libraries(${EXEC_NAME_0} ${INVARIANT_LIB_NAME} ${EXTERNAL_LIBS} ${LIB_NAME_1}) 
