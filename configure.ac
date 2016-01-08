dnl Process this file with autoconf to produce a configure script.

AC_INIT(elaboradar, [0.1], [palberoni@arpa.emr.it])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_MACRO_DIR([m4])

LT_INIT()
LT_LANG([C++])
AM_PROG_LIBTOOL

dnl To use subdirs
AC_PROG_MAKE_SET

AC_LANG([C++])

AC_ISC_POSIX
AC_PROG_CC
AC_PROG_CXX

dnl Needed by subdir-objects
AM_PROG_CC_C_O

dnl Use c++11
AX_CXX_COMPILE_STDCXX_11

dnl Give me warnings
AX_CXXFLAGS_WARN_ALL

dnl Mandatory librarires
AC_PATH_PROG(LOG4C_CONFIG, [log4c-config])
if test x"$LOG4C_CONFIG" = x; then
	AC_MSG_ERROR([*** missing Log4C library])
fi
LOG4C_CFLAGS=$($LOG4C_CONFIG --cflags)
LOG4C_LIBS=$($LOG4C_CONFIG --libs)
AC_SUBST(LOG4C_CFLAGS)
AC_SUBST(LOG4C_LIBS)

AC_PATH_PROG(GSL_CONFIG, [gsl-config])
if test x"$GSL_CONFIG" = x; then
	AC_MSG_ERROR([*** missing GSL (GNU Scientific Library)])
fi
GSL_CFLAGS=$($GSL_CONFIG --cflags)
GSL_LIBS=$($GSL_CONFIG --libs)
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_PATH_PROG(GDAL_CONFIG, [gdal-config])
if test x"$GDAL_CONFIG" = x; then
	AC_MSG_ERROR([*** missing GDAL (Geospatial Data Abstraction Library)])
fi
GDAL_CFLAGS=$($GDAL_CONFIG --cflags)
GDAL_LIBS=$($GDAL_CONFIG --libs)
AC_SUBST(GDAL_CFLAGS)
AC_SUBST(GDAL_LIBS)

PKG_CHECK_MODULES([HDF5], [hdf5], [have_hdf5=yes], [have_hdf5=no])
if test $have_hdf5 = yes
then
        PKG_CHECK_MODULES([HDF5], [hdf5])
else
	AC_CHECK_LIB([hdf5], [H5Fis_hdf5], [have_hdf5=yes], [have_hdf5=no])
	AC_CHECK_HEADER([hdf5.h], [have_hdf5=yes], [have_hdf5=no])
	if test x"$have_hdf5" = x"no"; then
		AC_MSG_ERROR("*** missing hdf5")
	else
		HDF5_CFLAGS=""
		HDF5_LIBS="-lhdf5"
		AC_SUBST(HDF5_CFLAGS)
		AC_SUBST(HDF5_LIBS)
	fi

        AC_MSG_WARN([compiling without libhdf5])
fi
PKG_CHECK_MODULES([RADARLIB], [radarlib])
PKG_CHECK_MODULES([EIGEN3], [eigen3])

AC_CONFIG_FILES([
Makefile
elaboradar/Makefile
src/Makefile
testdata/Makefile
esplosione/Makefile
elaboradar.pc
])
AC_OUTPUT
