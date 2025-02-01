cmake_minimum_required(VERSION 3.8)
project(Iterative-HFold)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -flto")

SET(SOURCE
  src/W_final.cpp 
  src/pseudo_loop.cpp 
  src/Iterative-HFold.cpp
  src/cmdline.cc
  src/Result.cc	
  src/s_energy_matrix.cpp	
  src/Hotspot.cc
  src/sparse_tree.cc
)

set(constraints_SOURCE
    src/ViennaRNA/constraints/constraints.c
    src/ViennaRNA/constraints/hard.c
    src/ViennaRNA/constraints/SHAPE.c
    src/ViennaRNA/constraints/soft.c
) 
set(datastructures_SOURCE
    src/ViennaRNA/datastructures/basic_datastructures.c
    src/ViennaRNA/datastructures/char_stream.c
    src/ViennaRNA/datastructures/lists.c
)
set(io_SOURCE
    src/ViennaRNA/io/file_formats.c
    src/ViennaRNA/io/io_utils.c
)
set(landscape_SOURCE
    src/ViennaRNA/landscape/move.c
)
set(loops_SOURCE
    src/ViennaRNA/loops/external_pf.c
    src/ViennaRNA/loops/external.c
    src/ViennaRNA/loops/hairpin.c
    src/ViennaRNA/loops/internal.c
    src/ViennaRNA/loops/multibranch.c
)
set(params_SOURCE
    src/ViennaRNA/params/default.c
    src/ViennaRNA/params/params.c
    src/ViennaRNA/params/io.c
)
set(utils_SOURCE
    src/ViennaRNA/utils/cpu.c
    src/ViennaRNA/utils/higher_order_functions.c
    src/ViennaRNA/utils/string_utils.c
    src/ViennaRNA/utils/structure_utils.c
    src/ViennaRNA/utils/utils.c
)
set(vienna_SOURCE
    ${constraints_SOURCE} ${datastructures_SOURCE}
    ${io_SOURCE} ${landscape_SOURCE}
    ${loops_SOURCE} ${params_SOURCE}
    ${utils_SOURCE}
    src/ViennaRNA/alphabet.c
    src/ViennaRNA/boltzmann_sampling.c
    src/ViennaRNA/centroid.c
    src/ViennaRNA/cofold.c
    src/ViennaRNA/commands.c
    src/ViennaRNA/dp_matrices.c
    src/ViennaRNA/equilibrium_probs.c
    src/ViennaRNA/eval.c
    src/ViennaRNA/fold_compound.c
    src/ViennaRNA/fold.c
    src/ViennaRNA/gquad.c
    src/ViennaRNA/grammar.c
    src/ViennaRNA/MEA.c
    src/ViennaRNA/mfe.c
    src/ViennaRNA/mm.c
    src/ViennaRNA/model.c
    src/ViennaRNA/part_func.c
    src/ViennaRNA/ribo.c
    src/ViennaRNA/sequence.c
    src/ViennaRNA/subopt.c
    src/ViennaRNA/unstructured_domains.c  
)

add_library(RNA STATIC ${vienna_SOURCE})
# configure_file(RNA.pc.in RNA.pc @ONLY)

target_include_directories(RNA PRIVATE .)

add_executable(Iterative-HFold ${SOURCE})

target_link_libraries(Iterative-HFold PRIVATE RNA)

include_directories(src)
