cmake_minimum_required(VERSION 3.25.1)

project(scalar2to2)

# Define the list of external libraries
set(EXTERNAL_LIBS
    pthread 
    looptools 
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
set(GENERIC_FLAGS "-Wall -Wextra -Wpedantic -Wno-deprecated-declarations -fPIC -fdiagnostics-color=always")
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

set(SOURCE_LIB_DIR src)
set(SOURCE_CLIB_DIR clib)

# Automatically gather all .cpp files in the src/ directory
file(GLOB LIB_CXX_SOURCES "${SOURCE_LIB_DIR}/*.cpp")

# Automatically gather all .c files in the src/ directory
file(GLOB LIB_C_SOURCES "${SOURCE_CLIB_DIR}/*.c")

# Automatically gather all .h files in the src/ directory
file(GLOB LIB_C_HEADERS "${SOURCE_CLIB_DIR}/*.h")

# Add the library and include all the .cpp files
add_library(${PROJECT_NAME} ${LIB_CXX_SOURCES} ${LIB_C_SOURCES})

# Specify output folder for the libraries
set_target_properties(${PROJECT_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Specify the include directory
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include  # C++ lib files
    ${CMAKE_CURRENT_SOURCE_DIR}/clib     # C lib files
)

# add dynamically linked libraries to the target library
target_link_libraries(${PROJECT_NAME} PRIVATE ${EXTERNAL_LIBS})
target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILE_LINK_FLAGS})

# Copying the headers of the C++ lib
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
# Copying the headers of the C lib
install(FILES ${LIB_C_HEADERS} DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
