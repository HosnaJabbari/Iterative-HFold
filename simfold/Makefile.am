## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS= subdir-objects

ACLOCAL_AMFLAGS = -I m4

AM_CXXFLAGS += -I$(srcdir)/include -I$(srcdir)/src/simfold -I$(srcdir)/src/common

lib_LTLIBRARIES=libsimfold.la

libsimfold_la_SOURCES = src/common/params.cpp src/common/common.cpp	\
  src/common/init.cpp src/common/timer.cpp				\
  src/simfold/s_specific_functions.cpp src/simfold/s_min_folding.cpp	\
  src/simfold/s_hairpin_loop.cpp src/simfold/s_stacked_pair.cpp		\
  src/simfold/s_energy_matrix.cpp src/simfold/s_internal_loop.cpp	\
  src/simfold/s_multi_loop.cpp src/simfold/s_multi_loop_sub.cpp		\
  src/simfold/s_sub_folding.cpp src/simfold/s_partition_function.cpp 

libsimfold_includedir = $(includedir)
libsimfold_include_HEADERS = src/simfold/s_min_folding.h		\
  src/simfold/s_multi_loop_sub.h					\
  src/simfold/s_partition_function_RESTRICTED.h				\
  src/simfold/s_partition_function_UNRESTRICTED.h			\
  src/simfold/s_specific_functions.h					\
  src/simfold/s_partition_function_complex.h				\
  src/simfold/s_stacked_pair.h src/simfold/s_hairpin_loop.h		\
  src/simfold/s_partition_function.h src/simfold/s_internal_loop.h	\
  src/simfold/s_multi_loop.h src/simfold/s_sub_folding.h		\
  src/simfold/s_energy_matrix.h src/common/timer.h			\
  src/common/common.h src/common/params.h include/constants.h		\
  include/externs.h include/globals.h include/init.h			\
  include/simfold.h include/structs.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = simfold.pc

########################################
## data

paramsdir = $(prefix)/share/simfold/params
params_DATA = params/turner-tstacki-rna.dh				\
  params/mathews-triloop-dna.dh params/mathews-dangle-dna.dh		\
  params/turner-tstacki-rna.dg params/turner-asint1x2-rna.dat		\
  params/parameters_DP09_chopped.txt params/mathews-tloop-dna.dat	\
  params/turner-special-hairpins-rna-with-T99.dat			\
  params/mathews-stack-dna.dh params/turner-tloop-rna.dg		\
  params/turner-sint2-rna.dat params/mathews-tstackh-dna.dh		\
  params/turner-tstackh-rna.dat params/mathews-miscloop-dna.dh		\
  params/parameters_DP09.txt params/turner-miscloop-rna.dh		\
  params/turner-stack-rna.dg params/turner-sint2-rna.dg			\
  params/turner-tloop-rna.dh params/turner-sint4-rna.dat		\
  params/turner-special-hairpins-rna.dat params/mathews-sint4-dna.dh	\
  params/turner_partypes_fm363.txt params/turner-triloop-rna.dh		\
  params/turner-miscloop-rna.dat params/mathews-tstacki-dna.dat		\
  params/turner-triloop-rna.dat						\
  params/CG_best_parameters_ISMB2007.txt params/turner-loop-rna.dh	\
  params/mathews-sint4-dna.dat params/mathews-miscloop-dna.dat		\
  params/mathews-sint2-dna.dat params/mathews-dangle-dna.dat		\
  params/mathews-tstacki-dna.dh						\
  params/turner_parameters_fm363_dangles0.txt				\
  params/mathews-sint2-dna.dh params/turner-tstacki-rna.dat		\
  params/mathews-stack-dna.dat params/turner-dangle-rna.dat		\
  params/turner-sint4-rna.dh params/turner-sint4-rna.dg			\
  params/mathews-tstackh-dna.dat params/turner-miscloop-rna.dg		\
  params/turner-sint2-rna.dh params/mathews-tloop-dna.dh		\
  params/multirnafold.conf params/mathews-loop-dna.dh			\
  params/turner-dangle-rna.dh params/turner-tloop-rna.dat		\
  params/turner-triloop-rna.dg params/CG_parameters_v1.1.txt		\
  params/turner-asint1x2-rna.dg params/turner-tstackh-rna.dg		\
  params/mathews-asint1x2-dna.dh params/turner-stack-rna.dat		\
  params/turner-loop-rna.dat params/mathews-triloop-dna.dat		\
  params/turner-tstackh-rna.dh params/turner-loop-rna.dg		\
  params/turner-stack-rna.dh params/mathews-asint1x2-dna.dat		\
  params/turner_parameters_fm363_constrdangles.txt			\
  params/mathews-loop-dna.dat params/turner-dangle-rna.dg		\
  params/turner-asint1x2-rna.dh


########################################
## binaries

# link binaries to libsimfold
LDADD = libsimfold.la

bin_PROGRAMS = simfold simfold_pf test_get_counts

simfold_SOURCES=simfold.cpp
simfold_pf_SOURCES=simfold_pf.cpp
test_get_counts_SOURCES=test_get_counts.cpp
