# Process this file with automake to produce Makefile.in

TESTFILES = \
  0termico.prev \
  201403010915h_ray.txt \
  201403010915_temp.Ngrad \
  stat_bloc \
  test_vpr \
  vpr_hmax_GAT \
  2014-03-01-09-15-00.itgat.PVOL.0.h5 \
  201403010915mat_bloc.bin \
  last_vpr_GAT \
  stat_elevchange \
  vpr_GAT \
  201403010915h_rayinf.txt \
  201403010915mat_el.bin \
  stat_anap \
  temperature.txt \
  vpr_heat_GAT

EXTRA_DIST = $(addsuffix .xz,$(TESTFILES))

%: %.xz
	xz -cd $^ > $@

BUILT_SOURCES = $(TESTFILES)

clean-local:
	rm -f $(TESTFILES)
