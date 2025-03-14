# Ian Wark July 6 2017

SET(SOURCE
  simfold.cpp
  src/common/params.cpp src/common/common.cpp	
  src/common/init.cpp 		
  src/simfold/s_specific_functions.cpp src/simfold/s_min_folding.cpp	
  src/simfold/s_hairpin_loop.cpp src/simfold/s_stacked_pair.cpp		
  src/simfold/s_energy_matrix.cpp src/simfold/s_internal_loop.cpp	
  src/simfold/s_multi_loop.cpp src/simfold/s_multi_loop_sub.cpp		
  src/simfold/s_sub_folding.cpp src/simfold/s_partition_function.cpp
  src/simfold/Hotspot.cc
  ../src/param_path.cc
)

SET(SOURCE_pf
  simfold_pf.cpp
  src/common/params.cpp src/common/common.cpp	
  src/common/init.cpp 		
  src/simfold/s_specific_functions.cpp src/simfold/s_min_folding.cpp	
  src/simfold/s_hairpin_loop.cpp src/simfold/s_stacked_pair.cpp		
  src/simfold/s_energy_matrix.cpp src/simfold/s_internal_loop.cpp	
  src/simfold/s_multi_loop.cpp src/simfold/s_multi_loop_sub.cpp		
  src/simfold/s_sub_folding.cpp src/simfold/s_partition_function.cpp
  src/simfold/Hotspot.cc
  ../src/param_path.cc
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# set location to make .a file
# set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# set location to make executable file
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# create library at LIBRARY_OUTPUT_PATH
add_library (simfold ${SOURCE})
add_executable(simfold_exec ${SOURCE})
add_executable(simfold_pf ${SOURCE_pf})
set_target_properties(simfold_exec PROPERTIES OUTPUT_NAME simfold)

# Enable compiler to find files for simfold library
target_include_directories (simfold PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (simfold PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories (simfold PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/common)
target_include_directories (simfold PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/simfold)

target_include_directories (simfold_exec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (simfold_exec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories (simfold_exec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/common)
target_include_directories (simfold_exec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/simfold)


target_include_directories (simfold_pf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (simfold_pf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories (simfold_pf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/common)
target_include_directories (simfold_pf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/simfold)

