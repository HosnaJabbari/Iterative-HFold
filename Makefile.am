ACLOCAL_AMFLAGS  = -I m4

bin_PROGRAMS= HFold HFold_pkonly HFold_iterative

AM_CFLAGS  = -I./simfold/include -I./simfold/src/common -I./simfold/src/simfold -I ./ -I/usr/include/boost -w -std=c++11

AM_LDFLAGS  = -L./ -lsimfold -lm -lpthread -L/usr/include/boost -lboost_regex

# Sources to compile HFold_iterative from
HFold_iterative_SOURCES = \
	W_final.h W_final.cpp \
	pseudo_loop.h pseudo_loop.cpp \
	V_final.h V_final.cpp \
	VM_final.h VM_final.cpp \
	h_common.h h_common.cpp \
	hfold_validation.h hfold_validation.cpp \
	hfold_iterative.h HFold_iterative.cpp

# Sources to compile HFold_iterative from
HFold_SOURCES = \
	W_final.h W_final.cpp \
	pseudo_loop.h pseudo_loop.cpp \
	V_final.h V_final.cpp \
	VM_final.h VM_final.cpp \
	h_common.h h_common.cpp \
	hfold_validation.h hfold_validation.cpp \
	hfold.h HFold.cpp

# Sources to compile HFold_iterative from
HFold_pkonly_SOURCES = \
	W_final.h W_final.cpp \
	pseudo_loop.h pseudo_loop.cpp \
	V_final.h V_final.cpp \
	VM_final.h VM_final.cpp \
	h_common.h h_common.cpp \
	hfold_validation.h hfold_validation.cpp \
	hfold_pkonly.h HFold_pkonly.cpp