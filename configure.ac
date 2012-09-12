AC_INIT(HiPLARM, version-0.1)

#AC_PROG_CC

AC_ARG_WITH([lapack],
    [AC_HELP_STRING([--with-lapack=LIBS],
        [LAPACK shared library @<:@-llapack@:>@])],
   [LIBS="${LIBS} $withval"] [lflag=1], [LIBS="${LIBS} -llapack"] [lflag=0]
)

if test $lflag -eq 0; then
	  AC_MSG_ERROR([The user should link a version of lapack preferably an optimised version such as ATLAS, OpenBLAS, Intel MKL or AMD ACML. Please see the Install guide for more details.])
fi

#################### PLASMA Libraries ####################

AC_ARG_WITH([plasma-lib],
    [AC_HELP_STRING([--with-plasma-lib=DIR],
        [full path to PLASMA shared libraries @<:@@:>@])],
    [LDFLAGS="${LDFLAGS} -L$withval/lib"] [pflag=1]
    CFLAGS="${CFLAGS} -I$withval/include -DHIPLAR_WITH_PLASMA"
    PLIBS="-lcoreblas -lquark -lplasma", [pflag=0]
)
################# CUDA Lib   ########################

AC_ARG_WITH([cuda-home],
	[AC_HELP_STRING([--with-cuda-home=DIR],
		[full path to CUDA installation @<:@/usr/local/cuda@:>@])],
	[CUDA_HOME="$withval"]
    LDFLAGS="${LDFLAGS} -L${CUDA_HOME}/lib64"
    CFLAGS="${CFLAGS} -I${CUDA_HOME}/include"
    CLIBS="-lcuda -lcudart -lcublas "
)



#################### MAGMA Libraries ####################

AC_ARG_WITH([magma-lib],
	[AC_HELP_STRING([--with-magma-lib=DIR],
		[full path to MAGMA shared libraries @<:@@:>@])],
	[LDFLAGS="${LDFLAGS} -L$withval/lib"] [mflag=1]
    CFLAGS="${CFLAGS} -I$withval/include -I$withval/testing -DHIPLAR_WITH_MAGMA"
    MLIBS="-lmagma -lmagmablas", [mflag=0]
)

LIBS="${LIBS} ${MLIBS} ${CLIBS} ${PLIBS}"

if test $pflag -eq 0; then
	  if test $mflag -eq 0; then
	  AC_MSG_ERROR([MAGMA or PLASMA not found. Please install PLASMA, MAGMA or both and follow the instructions. Please see the Install guide for more details.])
	  fi	  
fi


: ${R_HOME=`R RHOME`}
AC_ARG_WITH([r-home],
    [AC_HELP_STRING([--with-r-home=DIR],
        [full path to R installation @<:@@:>@])],
    [R_HOME=$withval]
)
AC_MSG_RESULT([setting R home directory to ${R_HOME}])
AC_CHECK_PROG([check_R], [R], [yes],, [${R_HOME}/bin])
if test -z "${check_R}"; then
    AC_MSG_ERROR([R not found!])
fi

AC_MSG_RESULT([setting R header files to ${R_INCLUDE_DIR}])
CC=`"${R_HOME}/bin/R" CMD config CC`
R_CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
CFLAGS="${CFLAGS} ${R_CPICFLAGS} -I${R_INCLUDE_DIR}"

AC_SUBST([CFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([LIBS])

AC_OUTPUT(src/Makevars)
