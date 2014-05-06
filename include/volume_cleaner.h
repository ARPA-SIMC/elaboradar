#ifndef ARCHIVIATORE_VOLUME_CLEANER_H
#define ARCHIVIATORE_VOLUME_CLEANER_H

#include <vector>
#include <cstdio>

namespace cumbac {

unsigned get_new_cell_num(unsigned orig_cell_num, unsigned max_range=0);

struct Beams
{
  unsigned char data_z[1024];
  unsigned char data_d[1024];
  unsigned char data_v[1024];
  unsigned char data_w[1024];
};

struct BeamCleaner
{
  unsigned char bin_wind_magic_number;
  unsigned min_segment_length, max_segment_length;

  BeamCleaner();

  void print_config(FILE* out) const;

  void clean_beams(Beams& b, unsigned beam_size, std::vector<bool>& corrected) const;
};

}

#endif
