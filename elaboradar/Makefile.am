# Process this file with automake to produce Makefile.in

EXTRA_DIST =

lib_LTLIBRARIES = libelaboradar.la

AM_CXXFLAGS=-std=c++11

AM_CPPFLAGS = -I$(top_srcdir) \
	      $(RADARLIB_CFLAGS) \
	      $(HDF5_CFLAGS) \
	      $(EIGEN3_CFLAGS) \
	      $(shell $(LOG4C_CONFIG) --cflags) \
	      $(shell $(GSL_CONFIG) --cflags) \
	      $(shell $(GDAL_CONFIG) --cflags) \
	      -I$(top_srcdir)/../libSP20/include \
              -I$(top_srcdir)/../util_Z-R/include \
		-fpermissive 

AM_LDFLAGS = $(RADARLIB_LIBS) $(EIGEN3_LIBS) $(HDF5_LIBS) -lhdf5_cpp \
	  $(shell $(LOG4C_CONFIG) --libs) \
	  $(shell $(GSL_CONFIG) --libs) \
	  $(shell $(GDAL_CONFIG) --libs) \
          -L$(top_srcdir)/../libSP20/lib/ -lSP20_utility \
          -L$(top_srcdir)/../util_Z-R/lib/ -lfunc_Z_R

elaboradarincludedir = $(includedir)/elaboradar

nobase_dist_elaboradarinclude_HEADERS = \
  logging.h \
  utils.h \
  matrix.h \
  image.h \
  volume.h \
  sp20.h \
  odim.h \
  elev_fin.h \
  algo/anaprop.h \
  algo/azimuth_resample.h \
  algo/cleaner.h \
  algo/dbz.h \
  algo/top.h \
  algo/utils.h

libelaboradar_la_SOURCES = \
  logging.cpp \
  utils.cpp \
  image.cpp \
  volume.cpp \
  sp20.cpp \
  odim.cpp \
  algo/anaprop.cpp \
  algo/azimuth_resample.cpp \
  algo/cleaner.cpp \
  algo/dbz.cpp

if DO_TESTS

check_PROGRAMS = testsuite

testsuite_SOURCES = \
  tests/test-azimutmap.cpp \
  tests/test-functions.cpp \
  tests/test-volume.cpp \
  tests/tut-main.cpp
testsuite_LDADD = libelaboradar.la $(WIBBLE_LIBS)

check-local:
	if [ -z "$(VERBOSE)" ]; then \
		$(DEBUGGER) ./testsuite $(TEST_ARGS) ; \
	else \
		LOG4C_RCPATH=../testdata/log4c_verbose $(DEBUGGER) ./testsuite $(TEST_ARGS) ; \
	fi

endif
