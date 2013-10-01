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

LDLIBS= $(LDLIBS_SP20) $(LDLIBS_ZR) $(LDLIBS_NR) -lm -llog4c

CC=cc 

CFLAGS= -I../include -Wall  -fstack-protector $(COMMON_DFLAGS_MEDIUM)   $(IFLAGS) -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_NR_C) -I$(INCLUDEDIR_Z_R)

CFLAGS_BB= -I../include -Wall  -fstack-protector $(COMMON_DFLAGS_SHORT)  $(IFLAGS) -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_NR_C) -I$(INCLUDEDIR_Z_R)

CFLAGS_BB_VPR=  -I../include -Wall -fstack-protector $(COMMON_DFLAGS_SHORT) $(ADDED_DFLAGS_VPR) $(IFLAGS) -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_NR_C) -I$(INCLUDEDIR_Z_R)

CFLAGS_BB_VPR_CLASS=  -I../include -Wall -fstack-protector $(COMMON_DFLAGS_SHORT) $(ADDED_DFLAGS_VPR_CLASS) $(IFLAGS) -I$(INCLUDEDIR_SP20) -I$(INCLUDEDIR_NR_C) -I$(INCLUDEDIR_Z_R)

OBJ_BB_VPR_CLASS= cum_bac_SP20_BB_VPR_CLASS.o setstat.o setwork.o func_Q3d.o Q_components.o 
OBJ_BB_VPR= cum_bac_SP20_BB_VPR.o setstat.o setwork.o func_Q3d.o Q_components.o 
OBJ_BB= cum_bac_SP20_BB.o setstat.o setwork.o func_Q3d.o Q_components.o 
OBJ= cum_bac_SP20.o 

all: install

cum_bac_SP20_m: $(OBJ)
	$(CC) -o $@ $^ $(LDLIBS) 

cum_bac_SP20.o:cum_bac_SP20.c $(INCLUDEDIR_SP20)/func_SP20read.h  
	$(CC) -c  $(CFLAGS) $< -o $@ 

cum_bac_SP20_BB: $(OBJ_BB)

cum_bac_SP20_BB.o:cum_bac_SP20_BB_VPR.c $(INCLUDEDIR_SP20)/func_SP20read.h ../include/MP_par.h ../include/vpr_par.h  
	$(CC) -c  $(CFLAGS_BB) $< -o $@

cum_bac_SP20_BB_VPR: $(OBJ_BB_VPR)

cum_bac_SP20_BB_VPR.o:cum_bac_SP20_BB_VPR.c $(INCLUDEDIR_SP20)/func_SP20read.h ../include/MP_par.h ../include/vpr_par.h 
	$(CC) -c  $(CFLAGS_BB_VPR) $< -o $@

cum_bac_SP20_BB_VPR_CLASS: $(OBJ_BB_VPR_CLASS)

cum_bac_SP20_BB_VPR_CLASS.o:cum_bac_SP20_BB_VPR.c $(INCLUDEDIR_SP20)/func_SP20read.h ../include/MP_par.h ../include/vpr_par.h ../include/par_class.h  
	$(CC) -c -g  $(CFLAGS_BB_VPR_CLASS) $< -o $@

setwork.o:setwork.c

setstat.o: setstat.c

func_Q3d.o:func_Q3d.c ../include/qual_par.h
	$(CC) -c $< $(CFLAGS) -o $@

Q_components.o:Q_components.c ../include/qual_par.h
	$(CC) -c $< $(CFLAGS) -o $@

clean: 
	$(RM) $(OBJ_BB_VPR) $(OBJ) $(OBJ_BB_VPR_CLASS) *.o $(BINDIR)/*

install: cum_bac_SP20_BB_VPR_CLASS cum_bac_SP20_BB_VPR cum_bac_SP20_BB cum_bac_SP20_m
	[ -d $(BINDIR) ] || mkdir -p $(BINDIR)
	 mv cum_bac_SP20_BB_VPR_CLASS $(BINDIR)
	 mv  cum_bac_SP20_BB_VPR  $(BINDIR)
	 mv  cum_bac_SP20_BB  $(BINDIR)
	 mv  cum_bac_SP20_m $(BINDIR)

%: %.o
	$(CC) -o $@ $^ $(LDLIBS)
