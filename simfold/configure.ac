dnl initial information about the project
AC_INIT([simfold],[1.0],[jabbari@ualberta.ca])

dnl check if the source folder is correct
AC_CONFIG_SRCDIR([./simfold.cpp])


dnl check for C++ preprocessor and compiler and the library compiler
AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_CC

dnl special dir for m4
AC_CONFIG_MACRO_DIR([m4])

dnl automake initialisation and check for minimal automake API version 1.9
AM_INIT_AUTOMAKE([1.9])

dnl init libtool
LT_INIT

dnl use the C++ compiler for the following checks
AC_LANG_CPLUSPLUS

dnl TURN ON C++11
AX_CXX_COMPILE_STDCXX_11

dnl Checks for header files.
AC_HEADER_STDC

dnl Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

dnl distribute additional compiler and linker flags
dnl --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])

dnl files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([simfold.pc])

dnl generate the final Makefile etc.
AC_OUTPUT([Makefile])
