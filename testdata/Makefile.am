# Process this file with automake to produce Makefile.in

TESTFILES = \
  DBP2_070120141530_GATTATICO \
  DBP2_060220140140_GATTATICO \
  201401231530_temp.Ngrad \
  201401071530mat_el.bin \
  201401071530mat_bloc.bin \
  201401071530h_ray.txt \
  201401071530h_rayinf.txt \
  temperature.txt \
  201402060140h_ray.txt \
  201402060140h_rayinf.txt \
  201402060140mat_bloc.bin \
  201402060140mat_el.bin \
  DBP2_060220140140_GATTATICO_mod \
  DBP2_020520141110_BOLOGNA \
  201405091240mat_bloc.bin \
  201405091240mat_el.bin \
  201405091240h_rayinf.txt \
  2014-05-09-12-40-00.itgat.PVOL.0.h5 \
  201405091240h_ray.txt \
  vpr/0termico.prev \
  vpr/2014-03-01-01-35-00.itgat.PVOL.0.h5 \
  vpr/201403010135h_rayinf.txt \
  vpr/201403010135h_ray.txt \
  vpr/201403010135mat_bloc.bin \
  vpr/201403010135mat_el.bin \
  vpr/last_vpr_GAT \
  vpr/temperature.txt \
  vpr/vpr_hmax_GAT

EXTRA_DIST = $(addsuffix .xz,$(TESTFILES))

%: %.xz
	xz -cd $^ > $@

BUILT_SOURCES = $(TESTFILES)

clean-local:
	rm -f $(TESTFILES)
