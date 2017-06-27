dnl initial information about the project
AC_INIT([HFold_iterative],[1.0],[jabbari@ualberta.ca])

dnl check if the source folder is correct
AC_CONFIG_SRCDIR([./HFold_iterative.cpp])

dnl check for C++ preprocessor and compiler and the library compiler
AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_CC

dnl special dir for m4
AC_CONFIG_MACRO_DIR([m4])

dnl automake initialisation and check for minimal automake API version 1.9
AM_INIT_AUTOMAKE([1.9])


dnl distribute additional compiler and linker flags
dnl --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])

dnl #################
dnl # Debug option
dnl 
AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
    AC_HELP_STRING(
        [--enable-debug],
        [enable debug data generation (def=no)]
    ),
    debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
    AC_DEFINE([DEBUG], [], [Define to 1 for turning on debugging])
    CXXFLAGS="$CXXFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG], [], [Define to 1 for turning off debugging])
fi

## check pkg-config
AC_CHECK_PROG([PKG_CONFIG_AVAILABLE],pkg-config,[true],[false])

dnl get path to simfold supplied as argument
AC_ARG_VAR([SIMFOLD_HOME],
    AC_HELP_STRING(
        [SIMFOLD_HOME],
        [Installation directory of the Simfold library]
      )
)

dnl ensure SIMFOLD_HOME argument is supplied
AS_IF([test "$SIMFOLD_HOME" = "" ],[
    AC_MSG_ERROR([Please set the environment variable SIMFOLD_HOME
to the installation directory of the simfold library.])
])

dnl if SIMFOLD_HOME is supplied, set PKG_CONFIG_PATH
AS_IF([test "$SIMFOLD_HOME" != "" ],[
    export PKG_CONFIG_PATH="$SIMFOLD_HOME/lib/pkgconfig:$PKG_CONFIG_PATH"
])

dnl configure simfold library via pkg-config
SIMFOLD_OK=false
PKG_CHECK_MODULES([SIMFOLD],[simfold >= 1.0],[SIMFOLD_OK=true], [SIMFOLD_PK=false])

dnl if unable to find simfold library, display error message
AS_IF([test "$SIMFOLD_OK" = "true"],[],[
    AC_MSG_ERROR([Cannot find the pkgconfig file of the simfold library.
Try setting the environment variable SIMFOLD_HOME
to the installation path of the simfold library.],
    [])
])

dnl pass SIMFOLD_HOME directory to HFHFold_iterativeold binary
AC_DEFINE_UNQUOTED(SIMFOLD_HOME,"$SIMFOLD_HOME")

dnl set flags to include -lsimfold
LIBS="$SIMFOLD_LIBS $LIBS"
AM_CXXFLAGS="$SIMFOLD_CFLAGS $AM_CXXFLAGS"

dnl set path to the simfold lib via rpath
AS_IF([test "$SIMFOLD_HOME" != "" ],[
    LDFLAGS="-Wl,-rpath,$SIMFOLD_HOME/lib"
])

dnl use the C++ compiler for the following checks
AC_LANG_CPLUSPLUS

dnl TURN ON C++11
AX_CXX_COMPILE_STDCXX_11

dnl Checks for header files.
AC_HEADER_STDC

dnl Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

dnl files to generate via autotools (.am or .in source files)
dnl AC_CONFIG_FILES([src/Makefile])

dnl generate the final Makefile etc.
AC_OUTPUT([Makefile])
