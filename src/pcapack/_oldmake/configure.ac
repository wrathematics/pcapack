## autoconf configure.ac && autoreconf -if

AC_PREREQ([2.69])
AC_INIT([pcapack], [0.1.0])
AC_CONFIG_SRCDIR([include/pcapack.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

