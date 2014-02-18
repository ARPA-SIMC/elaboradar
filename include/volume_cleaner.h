#ifndef ARCHIVIATORE_VOLUME_CLEANER_H
#define ARCHIVIATORE_VOLUME_CLEANER_H

#include <vector>

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
  int min_segment_length, max_segment_length;

  void clean_beams(Beams& b, unsigned beam_size, std::vector<bool>& corrected) const;
};

#endif
