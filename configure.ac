dnl Process this file with autoconf to produce a configure script.

AC_INIT(elaboradar, [0.6], [palberoni@arpa.emr.it])
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
dnl Debian still has a log4c without pkg-config support, see https://bugs.debian.org/834036
dnl PKG_CHECK_MODULES(LOG4C, [log4c])
AM_PATH_LOG4C
PKG_CHECK_MODULES(TCLAP, [tclap])
PKG_CHECK_MODULES(GSL, [gsl])
PKG_CHECK_MODULES(GDAL, [gdal], , [have_gdal=no])
if test x$have_gdal != xno
then
    AC_DEFINE([HAVE_GDAL], 1, [gdal is available])
fi
PKG_CHECK_MODULES(PROJ, [proj])

dnl Fedora 20 has no pkg-config for hdf5-devel
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
fi

PKG_CHECK_MODULES([RADARLIB], [radarlib])
PKG_CHECK_MODULES([EIGEN3], [eigen3])

AC_CONFIG_FILES([
Makefile
radarelab/Makefile
src/Makefile
testdata/Makefile
esplosione/Makefile
doc/Makefile
radarelab.pc
])
AC_OUTPUT
