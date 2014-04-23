/*
 * =====================================================================================
 *
 *       Filename:  volume_cleaner.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  18/02/2014 12:19:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <algorithm>
#include "volume_cleaner.h"

namespace cumbac {

using namespace std;

// TODO: toglierlo
#define MAX_DIM 512

unsigned get_new_cell_num(unsigned orig_cell_num, unsigned max_range)
{
    // Lunghezza che vogliamo
    unsigned rmax = min(orig_cell_num, (unsigned)MAX_DIM);
    if (max_range > 0)
      rmax = min(rmax, max_range);
    return rmax;
}


BeamCleaner::BeamCleaner()
  : min_segment_length(4), max_segment_length(40)
{
}

void BeamCleaner::print_config(FILE* out) const
{
    fprintf(out, "bin_wind_magic_number: %u, min_segment_length: %u, max_segment_length: %u\n",
      (unsigned)bin_wind_magic_number, min_segment_length, max_segment_length);
}

void BeamCleaner::clean_beams(Beams& b, unsigned beam_size, vector<bool>& corrected) const
{
/*----------------------------------------------------------------
INIZIO CODICE RIPULITURA 
-------------------------------------------------------------------*/
  int ibin;
  int start, end;
  int segment_length;
  bool in_a_segment = false;
  bool before, after;
  int ib, ia;
  for (ibin =0; ibin < beam_size; ibin++)
  {
    if (!in_a_segment) {
/* cerco la prima cella segmento da pulire*/
      //if (b.data_w[ibin] == 0 && b.data_v_w[ibin] == -125 )
      if (b.data_w[ibin] == 0 && b.data_v[ibin] == bin_wind_magic_number ){
        in_a_segment = true;
        start = ibin;
        after = false;
        before = false;
      }
    } else {
/* cerco la fine segmento da pulire*/
      //if (b.data_w[ibin] != 0 || b.data_v_w[ibin] != -125 || ibin == (beam_info.cell_num -1) )
      if (b.data_w[ibin] != 0 || b.data_v[ibin] != bin_wind_magic_number || ibin == (beam_size - 1) ){
        in_a_segment = false;
        end = ibin -1;
        if ( ibin == (beam_size -1) ) end=ibin;  // caso particolare per fine raggio
/* Fine trovata ora procedo alla pulizia eventuale */
        segment_length = end - start;
/* Cerco dati validi in Z prima del segmento */
        for (ib = ibin - 12; ib < ibin ; ib ++) 
          if (ib >= 0 && b.data_z[ib] >= 2 ) before = true;
/* Cerco dati validi in Z dopo il segmento */
        for (ia = ibin+1; ia <= ibin+12 ; ia ++) 
          if (ia < beam_size  && b.data_z[ia] >= 2 ) after = true;

        if ((segment_length >= min_segment_length && !before && !after) || 
             segment_length >= max_segment_length) {
/* qui pulisco */
   //         printf (" pulisco %d %d %d \n",segment_length, min_segment_length, max_segment_length);
          for (ib = start; ib <= end; ib++){
            b.data_z[ib] = 0;
            b.data_w[ib] = 0;
            *((signed char*)b.data_v + ib) = -128;
            corrected[ib] = true;
          }
        }
      }
    }
  }
}

}
