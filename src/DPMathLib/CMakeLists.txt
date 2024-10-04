cmake_minimum_required(VERSION 3.10)

project(DPMathLib)

# Define the list of external libraries
set(EXTERNAL_LIBS
    pthread 
    gsl 
    gslcblas
)
# Define the list of compile options (if needed)
set(COMPILE_LINK_FLAGS
    -lpthread
    -lgsl
    -lgslcblas
)

# Controlla se questo è il progetto principale
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    # Compilazione indipendente: usa i flag specifici del sottoprogetto

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
endif()

# Automatically gather all .cpp files in the src/ directory
file(GLOB LIB_CXX_SOURCES ${PROJECT_SOURCE_DIR}/*.cpp)

# Automatically gather all .c files in the src/ directory
file(GLOB LIB_C_SOURCES ${PROJECT_SOURCE_DIR}/*.c)

# Automatically gather all header files
file(GLOB LIB_HEADERS ${PROJECT_SOURCE_DIR}/*.h 
                      ${PROJECT_SOURCE_DIR}/*.hpp)

# Add the library and include all the .cpp files
add_library(${PROJECT_NAME} ${LIB_CXX_SOURCES} ${LIB_C_SOURCES})

# Specify output folder for the libraries
set_target_properties(${PROJECT_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Specify the include directory
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Apply the libraries and flags to the target
target_link_libraries(${PROJECT_NAME} PRIVATE ${EXTERNAL_LIBS})
target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILE_LINK_FLAGS})

# Installing header files
install(FILES ${LIB_HEADERS} DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
