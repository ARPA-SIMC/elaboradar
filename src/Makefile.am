# Process this file with automake to produce Makefile.in

EXTRA_DIST =

noinst_LTLIBRARIES = libtmp.la

AM_CXXFLAGS=-std=c++11

AM_CPPFLAGS = -I$(top_srcdir) \
              -I$(top_srcdir)/include \
	      $(RADARLIB_CFLAGS) \
	      $(HDF5_CFLAGS) \
	      $(EIGEN3_CFLAGS) \
	      $(LOG4C_CFLAGS) \
	      $(GSL_CFLAGS) \
	      $(GDAL_CFLAGS) \
	      -fpermissive

AM_LDFLAGS = $(RADARLIB_LIBS) $(EIGEN3_LIBS) $(HDF5_LIBS) -lhdf5_cpp $(LOG4C_LIBS) $(GSL_LIBS) $(GDAL_LIBS)


libtmp_la_SOURCES = \
  algo_steiner.cpp \
  algo_viz.cpp \
  config.cpp \
  cylindrical.cpp \
  interpola_vpr.cpp \
  interpola_vpr_gsl.cpp \
  interpola_vpr_nr.cpp \
  setwork.cpp \
  func_Q3d.cpp \
  Q_components.cpp \
  cart.cpp \
  cartproducts.cpp \
  site.cpp \
  assets.cpp \
  cum_bac.cpp \
  cum_bac_clparser.cpp
libtmp_la_LIBADD = ../elaboradar/libelaboradar.la

bin_PROGRAMS = elaboradar classificatore

elaboradar_SOURCES = \
  elaboradar.cpp
elaboradar_LDADD = ../elaboradar/libelaboradar.la libtmp.la

classificatore_SOURCES = \
  classificatore.cpp \
  classifier.cpp \
  melting_layer.cpp
classificatore_LDADD = ../elaboradar/libelaboradar.la libtmp.la
  

if DO_TESTS

check_PROGRAMS = ertest

ertest_SOURCES = \
  test-utils.cpp \
  tests/test-assets.cpp \
  tests/test-interpola-vpr.cpp \
  tests/test-cart.cpp \
  tests/test-processing.cpp \
  tests/test-read.cpp \
  tests/test-utils.cpp \
  tests/tut-main.cpp
ertest_LDADD = ../elaboradar/libelaboradar.la libtmp.la $(WIBBLE_LIBS)

check-local:
	if [ -z "$(VERBOSE)" ]; then \
		$(DEBUGGER) ./ertest $(TEST_ARGS) ; \
	else \
		LOG4C_RCPATH=../testdata/log4c_verbose $(DEBUGGER) ./ertest $(TEST_ARGS) ; \
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
