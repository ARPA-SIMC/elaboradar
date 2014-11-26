# Process this file with automake to produce Makefile.in

EXTRA_DIST =

lib_LTLIBRARIES = libelaboradar.la

AM_CXXFLAGS=-std=c++11

AM_CPPFLAGS = -I$(top_srcdir)/include \
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
  volume.h

libelaboradar_la_SOURCES = \
  volume.cpp

if DO_TESTS

check_PROGRAMS = testsuite

testsuite_SOURCES = \
  tests/tut-main.cpp
testsuite_LDADD = libelaboradar.la $(WIBBLE_LIBS)

check-local:
	if [ -z "$(VERBOSE)" ]; then \
		$(DEBUGGER) ./testsuite $(TEST_ARGS) ; \
	else \
		LOG4C_RCPATH=testdata/log4c_verbose $(DEBUGGER) ./testsuite $(TEST_ARGS) ; \
	fi

endif

## CFLAGS for all current code
#CFLAGS=-g --std=c++11 -I../include -Wall -O2 -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_Z_R) -I$(INCLUDEDIR_NR_C) \
#	   -DHAVE_INLINE \
#	   $(shell log4c-config --cflags) \
#	   $(shell pkg-config --cflags radarlib eigen3) \
#	   $(shell gsl-config --cflags) \
#	   $(shell gdal-config --cflags)
## Scommentare per disabilitare il range checking di GSL e guadagnare un po' in velocità
## CFLAGS += -DGSL_RANGE_CHECK_OFF
## Scommentare per disabilitare il range checking di Eigen e guadagnare un po' in velocità
## CFLAGS += -DEIGEN_NO_DEBUG
## Opzioni per debugging: -D_FORTIFY_SOURCE=2 -fstack-protector -Wstack-protector 
#
## CFLAGS for old code
#CFLAGS_LEGACY=$(CFLAGS) -DANAPROP -DTIME -DZ_LOWRIS -DMEDIUM
#
#CXXFLAGS=$(CFLAGS)
#
#OBJ=cum_bac_SP20.o
#
#
#all: install
#
#cum_bac_SP20_m: $(OBJ)
#	$(CXX) -o $@ $^ $(LDLIBS) 
#
#cum_bac_SP20.o:cum_bac_SP20.cpp $(INCLUDEDIR_SP20)/func_SP20read.h
#	$(CXX) -c  $(CFLAGS_LEGACY) $< -o $@ 
#
#elaboradar: elaboradar.o $(OBJ_COMMON)
#
#elaboradar.o: $(INCLUDEDIR_SP20)/func_SP20read.h ../include/MP_par.h ../include/vpr_par.h ../include/cum_bac_clparser.h
#
#setwork.o:setwork.cpp
#	$(CXX) -c $< $(CFLAGS) -o $@
#
#func_Q3d.o:func_Q3d.cpp ../include/qual_par.h
#	$(CXX) -c $< $(CFLAGS) -o $@
#
#Q_components.o:Q_components.cpp ../include/qual_par.h
#	$(CXX) -c $< $(CFLAGS) -o $@
#
#cum_bac_clparser.o:cum_bac_clparser.cpp ../include/cum_bac_clparser.h
#	$(CXX) -c $< $(CFLAGS) -o $@
#
#clean: 
#	$(RM) $(OBJ_COMMON) *.o $(BINDIR)/* tests/*.o cbtest
#	cd testdata && $(RM) $(TESTFILES)
#	cd esplosione && $(RM) $(TESTFILES_ESPLOSIONE)
#
#install: elaboradar  cum_bac_SP20_m
#	[ -d $(BINDIR) ] || mkdir -p $(BINDIR)
#	 mv  elaboradar  $(BINDIR)
#	 mv  cum_bac_SP20_m $(BINDIR)
#
## Compile test suite
#cbtest: $(patsubst %.cpp,%.o,$(wildcard tests/*.cpp)) $(OBJ_COMMON) test-utils.o
#	$(CXX) -o $@ $^ $(LDLIBS) `pkg-config --libs libwibble`
# 
#check: cbtest $(addprefix testdata/,$(TESTFILES)) $(addprefix esplosione/,$(TESTFILES_ESPLOSIONE))
#	if [ -z "$(VERBOSE)" ]; then \
#		$(DEBUGGER) ./cbtest $(TEST_ARGS) ; \
#	else \
#		LOG4C_RCPATH=testdata/log4c_verbose $(DEBUGGER) ./cbtest $(TEST_ARGS) ; \
#	fi
#
#%: %.o
#	$(CXX) -o $@ $^ $(LDLIBS)
#
