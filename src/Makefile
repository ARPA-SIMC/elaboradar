# se non è setata la libSP20_DIR vengono specificate qui

libSP20_DIR=../../libSP20/lib
libZR_DIR=../../util_Z-R/lib 
libNR_DIR=../../NR_C/lib

BINDIR=../bin

INCLUDEDIR_SP20=../../libSP20/include
INCLUDEDIR_NR_C=../../NR_C/include
INCLUDEDIR_Z_R=../../util_Z-R/include

LDLIBS_SP20= -L$(libSP20_DIR) -lSP20_utility
LDLIBS_ZR= -L$(libZR_DIR) -lfunc_Z_R
LDLIBS_NR=-L$(libNR_DIR) -lNR_C_ANSI

LDLIBS=$(LDLIBS_SP20) $(LDLIBS_ZR) $(LDLIBS_NR) -lm $(shell log4c-config --libs) $(shell pkg-config --libs radarlib) -lhdf5_cpp -lhdf5 $(shell gsl-config --libs) $(shell pkg-config eigen3 --libs)

CC=cc 
CXX=g++ 

# CFLAGS for all current code
CFLAGS=-g --std=c++11 -I../include -Wall -O2 -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_Z_R) -I$(INCLUDEDIR_NR_C) $(shell log4c-config --cflags) $(shell pkg-config --cflags radarlib) -DHAVE_INLINE $(shell gsl-config --cflags) $(shell pkg-config eigen3 --cflags)
# Scommentare per disabilitare il range checking di GSL e guadagnare un po' in velocità
# CFLAGS += -DGSL_RANGE_CHECK_OFF
# Scommentare per disabilitare il range checking di Eigen e guadagnare un po' in velocità
# CFLAGS += -DEIGEN_NO_DEBUG
# Opzioni per debugging: -D_FORTIFY_SOURCE=2 -fstack-protector -Wstack-protector 

# CFLAGS for old code
CFLAGS_LEGACY=$(CFLAGS) -DANAPROP -DTIME -DZ_LOWRIS -DMEDIUM

CXXFLAGS=$(CFLAGS)

OBJ_VOLUME=volume.o volume_cleaner.o volume_loader.o volume_sp20.o volume_odim.o

OBJ_COMMON=$(OBJ_VOLUME) cylindrical.o steiner.o viz.o cum_bac.o interpola_vpr.o interpola_vpr_gsl.o interpola_vpr_nr.o  setwork.o func_Q3d.o Q_components.o utils.o site.o assets.o logging.o cum_bac_clparser.o

OBJ=cum_bac_SP20.o

TESTFILES = DBP2_070120141530_GATTATICO DBP2_060220140140_GATTATICO 201401231530_temp.Ngrad 201401071530mat_el.bin 201401071530mat_bloc.bin 201401071530h_ray.txt 201401071530h_rayinf.txt temperature.txt 201402060140h_ray.txt 201402060140h_rayinf.txt 201402060140mat_bloc.bin 201402060140mat_el.bin DBP2_060220140140_GATTATICO_mod DBP2_020520141110_BOLOGNA 201405091240mat_bloc.bin 201405091240mat_el.bin 201405091240h_rayinf.txt 2014-05-09-12-40-00.itgat.PVOL.0.h5 201405091240h_ray.txt


all: install

cum_bac_SP20_m: $(OBJ)
	$(CXX) -o $@ $^ $(LDLIBS) 

cum_bac_SP20.o:cum_bac_SP20.cpp $(INCLUDEDIR_SP20)/func_SP20read.h  
	$(CXX) -c  $(CFLAGS_LEGACY) $< -o $@ 

elaboradar: elaboradar.o $(OBJ_COMMON)

elaboradar.o: $(INCLUDEDIR_SP20)/func_SP20read.h ../include/MP_par.h ../include/vpr_par.h ../include/cum_bac_clparser.h

setwork.o:setwork.cpp
	$(CXX) -c $< $(CFLAGS) -o $@

func_Q3d.o:func_Q3d.cpp ../include/qual_par.h
	$(CXX) -c $< $(CFLAGS) -o $@

Q_components.o:Q_components.cpp ../include/qual_par.h
	$(CXX) -c $< $(CFLAGS) -o $@

cum_bac_clparser.o:cum_bac_clparser.cpp ../include/cum_bac_clparser.h
	$(CXX) -c $< $(CFLAGS) -o $@

clean: 
	$(RM) $(OBJ_COMMON) *.o $(BINDIR)/* tests/*.o cbtest
	cd testdata && $(RM) $(TESTFILES)

install: elaboradar  cum_bac_SP20_m
	[ -d $(BINDIR) ] || mkdir -p $(BINDIR)
	 mv  elaboradar  $(BINDIR)
	 mv  cum_bac_SP20_m $(BINDIR)

# Compile test suite
cbtest: $(patsubst %.cpp,%.o,$(wildcard tests/*.cpp)) $(OBJ_COMMON) test-utils.o
	$(CXX) -o $@ $^ $(LDLIBS) `pkg-config --libs libwibble`
 
check: cbtest $(addprefix testdata/,$(TESTFILES))
	if [ -z "$(VERBOSE)" ]; then \
		$(DEBUGGER) ./cbtest $(TEST_ARGS) ; \
	else \
		LOG4C_RCPATH=testdata/log4c_verbose $(DEBUGGER) ./cbtest $(TEST_ARGS) ; \
	fi

%: %.o
	$(CXX) -o $@ $^ $(LDLIBS)

%: %.xz
	xz -cd $^ > $@
