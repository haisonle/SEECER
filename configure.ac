AC_INIT([seecer], [1.0], [hple@cs.cmu.edu,mschulz@cs.cmu.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
	src/Makefile
])

AC_OPENMP
AC_CHECK_LIB([gomp],[GOMP_parallel_start])

AC_CHECK_LIB([rt],[clock_gettime])
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AC_ARG_WITH([seqan-include-path],
  [AS_HELP_STRING([--seqan-include-path],
    [location of the SeqAn headers, default to ../../seqan-1.3.1/])],
  [SEQAN_CFLAGS="-I$withval"],
  [SEQAN_CFLAGS='-I../../seqan-1.3.1/'])
AC_SUBST(SEQAN_CFLAGS)

AC_OUTPUT
