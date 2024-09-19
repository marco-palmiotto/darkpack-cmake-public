cmake_minimum_required(VERSION 3.10)

project(test)

set(INVARIANT_LIB_NAME MathLib)

set(LIB_NAME_1 FinalLib)
set(EXEC_NAME_0 main_test.x)

# Using the installed path for the invariant library's header
set(INVARIANT_LIB_HEADER_PATH ${CMAKE_BINARY_DIR}/include/${INVARIANT_LIB_NAME})  

add_executable(${EXEC_NAME_0} main.cpp)

# for each executable, we perform some operations:
# 1. defining new header location
target_include_directories(${EXEC_NAME_0} PUBLIC
    ${INVARIANT_LIB_HEADER_PATH}
    ${CMAKE_BINARY_DIR}/include
)

# 2. manual link to the path of the libraries
target_link_directories(${EXEC_NAME_0} PRIVATE ${CMAKE_BINARY_DIR}/lib/)
# 3. manual link of the library
target_link_libraries(${EXEC_NAME_0} ${INVARIANT_LIB_NAME} ${LIB_NAME_1})
