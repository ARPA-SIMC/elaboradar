# se non è setata la libSP20_DIR vengono specificate qui

libSP20_DIR=../../libSP20/lib
libZR_DIR=../../util_Z-R/lib 
libNR_DIR=../../NR_C/lib

BINDIR=../bin

INCLUDEDIR_SP20=../../libSP20/include
INCLUDEDIR_NR_C=../../NR_C/include
INCLUDEDIR_Z_R=../../util_Z-R/include

COMMON_DFLAGS_MEDIUM= -D_FORTIFY_SOURCE=2 -DANAPROP -DTIME -DZ_LOWRIS  -DMEDIUM 
COMMON_DFLAGS_SHORT= -D_FORTIFY_SOURCE=2 -DANAPROP -DTIME  -DZ_LOWRIS  -DSHORT -DSTATIC  -DBEAMBLOCKING -DQUALITY 
ifeq ($(env),devel)
COMMON_DFLAGS_SHORT= -D_FORTIFY_SOURCE=2 -DANAPROP -DTIME  -DZ_LOWRIS  -DSHORT -DSTATIC  -DBEAMBLOCKING -DQUALITY -DSTAMPE_EXTRA
endif

ADDED_DFLAGS_VPR= -DVPR
ADDED_DFLAGS_VPR_CLASS= -DVPR -DCLASS

LDLIBS_SP20= -L$(libSP20_DIR) -lSP20_utility
LDLIBS_ZR= -L$(libZR_DIR) -lfunc_Z_R
LDLIBS_NR=-L$(libNR_DIR) -lNR_C_ANSI

LDLIBS= $(LDLIBS_SP20) $(LDLIBS_ZR) $(LDLIBS_NR) -lm $(shell log4c-config --libs)

CC=cc 
CXX=g++ 

CFLAGS=-g -I../include -Wall  -fstack-protector $(COMMON_DFLAGS_MEDIUM)   $(IFLAGS) -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_NR_C) -I$(INCLUDEDIR_Z_R) $(shell log4c-config --cflags)
CXXFLAGS=$(CFLAGS)

CFLAGS_BB=-g -I../include -Wall  -fstack-protector $(COMMON_DFLAGS_SHORT)  $(IFLAGS) -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_NR_C) -I$(INCLUDEDIR_Z_R)

CFLAGS_BB_VPR=-g -I../include -Wall -fstack-protector $(COMMON_DFLAGS_SHORT) $(ADDED_DFLAGS_VPR) $(IFLAGS) -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_NR_C) -I$(INCLUDEDIR_Z_R)

CFLAGS_BB_VPR_CLASS=-g -I../include -Wall -fstack-protector $(COMMON_DFLAGS_SHORT) $(ADDED_DFLAGS_VPR_CLASS) $(IFLAGS) -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_NR_C) -I$(INCLUDEDIR_Z_R)

OBJ_COMMON=setstat.o setwork.o func_Q3d.o Q_components.o utils.o logging.o
OBJ_BB_VPR_CLASS=cum_bac.o cum_bac_SP20_BB_VPR_CLASS.o $(OBJ_COMMON)
OBJ_BB_VPR=cum_bac.o cum_bac_SP20_BB_VPR.o $(OBJ_COMMON)
OBJ_BB=cum_bac.o cum_bac_SP20_BB.o $(OBJ_COMMON)
OBJ= cum_bac_SP20.o 

all: install

cum_bac_SP20_m: $(OBJ)
	$(CXX) -o $@ $^ $(LDLIBS) 

cum_bac_SP20.o:cum_bac_SP20.cpp $(INCLUDEDIR_SP20)/func_SP20read.h  
	$(CXX) -c  $(CFLAGS) $< -o $@ 

cum_bac_SP20_BB: $(OBJ_BB)

cum_bac_SP20_BB.o:cum_bac_SP20_BB_VPR.cpp $(INCLUDEDIR_SP20)/func_SP20read.h ../include/MP_par.h ../include/vpr_par.h  
	$(CXX) -c  $(CFLAGS_BB) $< -o $@

cum_bac_SP20_BB_VPR: $(OBJ_BB_VPR)

cum_bac_SP20_BB_VPR.o:cum_bac_SP20_BB_VPR.cpp $(INCLUDEDIR_SP20)/func_SP20read.h ../include/MP_par.h ../include/vpr_par.h 
	$(CXX) -c  $(CFLAGS_BB_VPR) $< -o $@

cum_bac_SP20_BB_VPR_CLASS: $(OBJ_BB_VPR_CLASS)

cum_bac_SP20_BB_VPR_CLASS.o:cum_bac_SP20_BB_VPR.cpp $(INCLUDEDIR_SP20)/func_SP20read.h ../include/MP_par.h ../include/vpr_par.h ../include/par_class.h  
	$(CXX) -c  $(CFLAGS_BB_VPR_CLASS) $< -o $@

setwork.o:setwork.cpp
	$(CXX) -c $< $(CFLAGS) -o $@

setstat.o: setstat.cpp
	$(CXX) -c $< $(CFLAGS) -o $@

func_Q3d.o:func_Q3d.cpp ../include/qual_par.h
	$(CXX) -c $< $(CFLAGS) -o $@

Q_components.o:Q_components.cpp ../include/qual_par.h
	$(CXX) -c $< $(CFLAGS) -o $@

clean: 
	$(RM) $(OBJ_BB_VPR) $(OBJ) $(OBJ_BB_VPR_CLASS) *.o $(BINDIR)/*

install: cum_bac_SP20_BB_VPR_CLASS cum_bac_SP20_BB_VPR cum_bac_SP20_BB cum_bac_SP20_m
	[ -d $(BINDIR) ] || mkdir -p $(BINDIR)
	 mv cum_bac_SP20_BB_VPR_CLASS $(BINDIR)
	 mv  cum_bac_SP20_BB_VPR  $(BINDIR)
	 mv  cum_bac_SP20_BB  $(BINDIR)
	 mv  cum_bac_SP20_m $(BINDIR)

# Compile test suite
cbtest: $(patsubst %.cpp,%.o,$(wildcard tests/*.cpp)) cum_bac.o $(OBJ_COMMON)
	$(CXX) -o $@ $^ $(LDLIBS) `pkg-config --libs libwibble`
 
check: cbtest testdata/DBP2_070120141530_GATTATICO
	./cbtest


%: %.o
	$(CXX) -o $@ $^ $(LDLIBS)

%: %.xz
	xz -cd $^ > $@
