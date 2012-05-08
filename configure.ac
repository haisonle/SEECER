AC_INIT([seecer], [1.0], [hple@cs.cmu.edu,mschulz@cs.cmu.edu])
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects  foreign])
AC_PREREQ(2.62)

: ${CXXFLAGS=""}
CFLAGS="$CFLAGS -I/opt/local/include"
CXXFLAGS="$CXXFLAGS -I/opt/local/include"
LDFLAGS="$LDFLAGS -L/opt/local/lib"
AC_PROG_CXX

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
	Makefile
])

AC_OPENMP
AC_CHECK_LIB([gomp],[GOMP_parallel_start])

AC_CHECK_LIB([rt],[clock_gettime])
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[main], [], [AC_ERROR([Missing library lgslcblas])])
AC_CHECK_LIB([gsl],[main], [], [AC_ERROR([Missing library lgsl])])
AC_ARG_WITH([seqan-include-path],
  [AS_HELP_STRING([--seqan-include-path],
    [location of the SeqAn headers, default to seqan-cxx/])],
  [SEQAN_CFLAGS="-I$withval"],
  [SEQAN_CFLAGS='-Iseqan-cxx/'])
AC_SUBST(SEQAN_CFLAGS)

AC_OUTPUT
