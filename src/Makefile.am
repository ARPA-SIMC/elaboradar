# Process this file with automake to produce Makefile.in

EXTRA_DIST =

lib_LTLIBRARIES = libelaboradar.la

AM_CPPFLAGS = -I$(top_srcdir)/include \
	      $(RADARLIB_CFLAGS) \
	      $(HDF5_CFLAGS) \
	      $(EIGEN3_CFLAGS) \
	      $(shell $(LOG4C_CONFIG) --cflags) \
	      $(shell $(GSL_CONFIG) --cflags) \
	      $(shell $(GDAL_CONFIG) --cflags) \
	      -I$(top_srcdir)/../libSP20/include \
              -I$(top_srcdir)/../util_Z-R/include

AM_LDFLAGS = $(RADARLIB_LIBS) $(EIGEN3_LIBS) $(HDF5_LIBS) -lhdf5_cpp \
	  $(shell $(LOG4C_CONFIG) --libs) \
	  $(shell $(GSL_CONFIG) --libs) \
	  $(shell $(GDAL_CONFIG) --libs) \
          -L$(top_srcdir)/../libSP20/lib/ -lSP20_utility \
          -L$(top_srcdir)/../util_Z-R/lib/ -lfunc_Z_R


libelaboradar_la_SOURCES = \
  volume.cpp \
  volume_azimuthmap.cpp \
  volume_loader.cpp \
  volume_sp20.cpp \
  volume_odim.cpp \
  volume_resample.cpp \
  algo_cleaner.cpp \
  algo_steiner.cpp \
  algo_viz.cpp \
  algo_anaprop.cpp \
  algo_dbz.cpp \
  logging.cpp \
  config.cpp \
  utils.cpp \
  cylindrical.cpp \
  interpola_vpr.cpp \
  interpola_vpr_gsl.cpp \
  interpola_vpr_nr.cpp \
  setwork.cpp \
  func_Q3d.cpp \
  Q_components.cpp \
  site.cpp \
  assets.cpp \
  cum_bac.cpp \
  cum_bac_clparser.cpp \
  image.cpp

bin_PROGRAMS = elaboradar

elaboradar_SOURCES = \
  elaboradar.cpp
elaboradar_LDADD = libelaboradar.la

if DO_TESTS

check_PROGRAMS = ertest

ertest_SOURCES = \
  tests/test-assets.cpp \
  tests/test-azimutmap.cpp \
  tests/test-functions.cpp \
  tests/test-interpola-vpr.cpp \
  tests/test-processing.cpp \
  tests/test-read.cpp \
  tests/test-utils.cpp \
  tests/test-volume.cpp \
  tests/tut-main.cpp
ertest_LDADD = libelaboradar.la $(WIBBLE_LIBS)

check-local:
	if [ -z "$(VERBOSE)" ]; then \
		$(DEBUGGER) ./ertest $(TEST_ARGS) ; \
	else \
		LOG4C_RCPATH=testdata/log4c_verbose $(DEBUGGER) ./ertest $(TEST_ARGS) ; \
	fi

TESTFILES = \
  testdata/DBP2_070120141530_GATTATICO \
  testdata/DBP2_060220140140_GATTATICO \
  testdata/201401231530_temp.Ngrad \
  testdata/201401071530mat_el.bin \
  testdata/201401071530mat_bloc.bin \
  testdata/201401071530h_ray.txt \
  testdata/201401071530h_rayinf.txt \
  testdata/temperature.txt \
  testdata/201402060140h_ray.txt \
  testdata/201402060140h_rayinf.txt \
  testdata/201402060140mat_bloc.bin \
  testdata/201402060140mat_el.bin \
  testdata/DBP2_060220140140_GATTATICO_mod \
  testdata/DBP2_020520141110_BOLOGNA \
  testdata/201405091240mat_bloc.bin \
  testdata/201405091240mat_el.bin \
  testdata/201405091240h_rayinf.txt \
  testdata/2014-05-09-12-40-00.itgat.PVOL.0.h5 \
  testdata/201405091240h_ray.txt

TESTFILES_ESPLOSIONE = \
  esplosione/0termico.prev \
  esplosione/201403010915h_ray.txt \
  esplosione/201403010915_temp.Ngrad \
  esplosione/stat_bloc \
  esplosione/test_vpr \
  esplosione/vpr_hmax_GAT \
  esplosione/2014-03-01-09-15-00.itgat.PVOL.0.h5 \
  esplosione/201403010915mat_bloc.bin \
  esplosione/last_vpr_GAT \
  esplosione/stat_elevchange \
  esplosione/vpr_GAT \
  esplosione/201403010915h_rayinf.txt \
  esplosione/201403010915mat_el.bin \
  esplosione/stat_anap \
  esplosione/temperature.txt \
  esplosione/vpr_heat_GAT

EXTRA_DIST += $(addsuffix .xz,$(TESTFILES)) $(addsuffix .xz,$(TESTFILES_ESPLOSIONE))

%: %.xz
	xz -cd $^ > $@

BUILT_SOURCES = $(TESTFILES) $(TESTFILES_ESPLOSIONE)

clean-local:
	rm -f $(TESTFILES)
	rm -f $(TESTFILES_ESPLOSIONE)

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
