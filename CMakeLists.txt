# Ian Wark July 6 2017
# Edited by Connie He July 7 2022
cmake_minimum_required(VERSION 3.1)
project(Iterative-HFold)

# enable c++11
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -std=c++17 \ -O3")

# create simfold library
add_subdirectory(simfold)

set(Iterative-HFold_SOURCE
    W_final.cpp 
    pseudo_loop.cpp 
    V_final.cpp 
    VM_final.cpp 
    h_common.cpp 
    hfold_validation.cpp 
    Iterative-HFold.cpp
    cmdline.cc
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# set path for executable to be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executables
add_executable(Iterative-HFold ${Iterative-HFold_SOURCE})

# link to simfold
target_link_libraries(Iterative-HFold LINK_PUBLIC simfold)

