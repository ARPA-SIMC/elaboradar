# Process this file with automake to produce Makefile.in

EXTRA_DIST =

lib_LTLIBRARIES = libradarelab.la

AM_CPPFLAGS = -I$(top_srcdir) \
	      $(RADARLIB_CFLAGS) \
	      $(HDF5_CFLAGS) \
	      $(EIGEN3_CFLAGS) \
	      $(LOG4C_CFLAGS) \
	      $(GSL_CFLAGS) \
	      $(GDAL_CFLAGS)

AM_LDFLAGS = $(RADARLIB_LIBS) $(EIGEN3_LIBS) $(HDF5_LIBS) -lhdf5_cpp $(LOG4C_LIBS) $(GSL_LIBS) $(GDAL_LIBS)

radarelabincludedir = $(includedir)/radarelab

nobase_dist_radarelabinclude_HEADERS = \
  logging.h \
  utils.h \
  utils/string.h \
  utils/sys.h \
  utils/subprocess.h \
  utils/term.h \
  utils/testrunner.h \
  utils/tests.h \
  cylindrical.h \
  matrix.h \
  image.h \
  interpola_vpr.h \
  loader.h \
  par_class.h \
  RadarSite.h \
  volume.h \
  vpr_par.h \
  sp20.h \
  odim.h \
  elev_fin.h \
  cart.h \
  algo/anaprop.h \
  algo/azimuth_resample.h \
  algo/cleaner.h \
  algo/dbz.h \
  algo/elabora_volume.h \
  algo/statistics.h \
  algo/steiner.h \
  algo/top.h \
  algo/viz.h \
  algo/vpr.h

libradarelab_la_SOURCES = \
  logging.cpp \
  utils.cpp \
  utils/string.cc \
  utils/sys.cc \
  utils/subprocess.cc \
  utils/term.cc \
  utils/testrunner.cc \
  utils/tests.cc \
  cylindrical.cpp \
  image.cpp \
  interpola_vpr.cpp \
  interpola_vpr_gsl.cpp \
  interpola_vpr_nr.cpp \
  volume.cpp \
  cart.cpp \
  sp20.cpp \
  sp20/SP20read.cpp \
  sp20/func_SP20toMDB.cpp \
  odim.cpp \
  algo/anaprop.cpp \
  algo/azimuth_resample.cpp \
  algo/cleaner.cpp \
  algo/dbz.cpp \
  algo/elabora_volume.cpp \
  algo/steiner.cpp \
  algo/viz.cpp \
  algo/vpr.cpp

check_PROGRAMS = testsuite

testsuite_SOURCES = \
  tests/test-azimutmap.cpp \
  tests/test-functions.cpp \
  tests/test-volume.cpp \
  tests/tests-main.cc
testsuite_LDADD = libradarelab.la $(WIBBLE_LIBS)

check-local:
	if [ -z "$(VERBOSE)" ]; then \
		libtool --mode=execute $(DEBUGGER) ./testsuite $(TEST_ARGS) ; \
	else \
		LOG4C_RCPATH=../testdata/log4c_verbose libtool --mode=execute $(DEBUGGER) ./testsuite $(TEST_ARGS) ; \
	fi
