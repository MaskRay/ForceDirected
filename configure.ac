AC_INIT([force-directed], [0.1], [emacs.ray@gmail.com])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

AC_LANG([C++])
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

AC_ARG_WITH([float], AC_HELP_STRING([--with-float], [Use float]))

AS_CASE(["$with_float"], [yes], [AC_DEFINE([USE_FLOAT], [], [Use float])])

AC_OUTPUT
