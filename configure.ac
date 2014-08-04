## autoconf configure.ac && autoreconf -if

AC_INIT(DESCRIPTION)

AC_CONFIG_SRCDIR([src/ffpca.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Obtain from R.
R_SCMD="${R_HOME}/bin/R CMD config"
FLIBS=`${R_SCMD} FLIBS`
BLAS_LIBS=`${R_SCMD} BLAS_LIBS`
LAPACK_LIBS=`${R_SCMD} LAPACK_LIBS`


# Check tools.
AC_CHECK_PROG([ls_ok], [ls], [yes], [no], ,)
AC_CHECK_PROG([echo_ok], [echo], [yes], [no], ,)
AC_CHECK_PROG([sed_ok], [sed], [yes], [no], ,)
AC_CHECK_PROG([grep_ok], [grep], [yes], [no], ,)


# Check Rlapack.
if test "X$echo_ok" = "Xyes" -a "X$grep_ok" = "Xyes"; then
  Rlapack_ok=`echo ${LAPACK_LIBS} | grep "\-lRlapack"`
  if test "X$Rlapack_ok" = "X"; then
    dnl Rlapack is not found.
    Rlapack_ok="NA"
  else
    echo "Found Rlpack from ${Rlapack_ok}"
  fi
else
  dnl echo and grep are not found.
  Rlapack_ok="NA"
fi

if test "X$Rlapack_ok" = "XNA"; then
  # An external LAPACK is used in run time.
  R_OBJS=""
else
  dnl Rlapack is used.

  ### XERBLA_ is not available.
  # AC_LANG(Fortran)
  # AC_SEARCH_LIBS([zgeev], [Rlapack], [zgeev_ok=yes], [zgeev_ok=no],
  # [$LAPACK_LIBS $BLAS_LIBS $FLIBS])
  # AC_CHECK_FUNC([zgeev], [zgeev_ok=yes], [zgeev_ok=no])

  ### F77_CALL() is not available.
  # R_INCLUDE="${R_HOME}/include/"
  # CPPFLAGS="-I${R_INCLUDE}"
  # AC_CHECK_HEADER([R.h])
  # AC_CHECK_HEADER([Rinternals.h])
  # AC_CHECK_HEADER([R_ext/Lapack.h])
  # AC_CHECK_FUNC([zgeev], [zgeev_ok=yes], [zgeev_ok=no])

  if test "X$ls_ok" = "Xyes" -a "X$sed_ok" = "Xyes" -a "X$grep_ok" = "Xyes"; then
    ### Check R_ext/Lapack.h
    LAPACK_HEADER="${R_HOME}/include/R_ext/Lapack.h"

    funcs=`ls src/pcapack/src/*.f | sed -e "s/.*\/\(.*\)\.f/\\1/"`
    R_OBJS=""
    for func in ${funcs}
    do
      func_ok=`grep "F77_NAME($func)" ${LAPACK_HEADER}`
      echo "Found $func from ${func_ok}"
      if test "X${func_ok}" = "X"; then
        R_OBJS="${R_OBJS} pcapack/$func.o"
      fi
    done
  else
    ### Use default if all above fail.
    R_OBJS="\$(R_OBJS)"
  fi
fi

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_SUBST(QZ_OBJS)
AC_OUTPUT(src/Makevars)

