# Ian Wark July 6 2017

cmake_minimum_required(VERSION 2.8)
project(HFold)

# only continue if we know we already have boost
# create simfold library
add_subdirectory(simfold)

set(HFold_SOURCE
    W_final.h W_final.cpp 
    pseudo_loop.h pseudo_loop.cpp 
    V_final.h V_final.cpp 
    VM_final.h VM_final.cpp
    h_common.h h_common.cpp 
    hfold_validation.h hfold_validation.cpp 
    hfold.h HFold.cpp
)

set(HFold_iterative_SOURCE
    W_final.h W_final.cpp 
    pseudo_loop.h pseudo_loop.cpp 
    V_final.h V_final.cpp 
    VM_final.h VM_final.cpp 
    h_common.h h_common.cpp 
    hfold_validation.h hfold_validation.cpp 
    hfold_iterative.h HFold_iterative.cpp
)

set(HFold_pkonly_SOURCE
    W_final.h W_final.cpp 
    pseudo_loop.h pseudo_loop.cpp 
    V_final.h V_final.cpp 
    VM_final.h VM_final.cpp 
    h_common.h h_common.cpp 
    hfold_validation.h hfold_validation.cpp 
    hfold_pkonly.h HFold_pkonly.cpp
)

# enable c++11
add_compile_options(-std=c++11)

# get pthread library
find_package (Threads)

# set path for executable to be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executables
add_executable(HFold ${HFold_SOURCE})
add_executable(HFold_iterative ${HFold_iterative_SOURCE})
add_executable(HFold_pkonly ${HFold_pkonly_SOURCE})

# link to simfold
target_link_libraries(HFold LINK_PUBLIC simfold)
target_link_libraries(HFold_iterative LINK_PUBLIC simfold)
target_link_libraries(HFold_pkonly LINK_PUBLIC simfold)

# link to pthread
target_link_libraries(HFold LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(HFold_iterative LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(HFold_pkonly LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT})