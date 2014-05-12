#ifndef ARCHIVIATORE_VOLUME_CLEANER_H
#define ARCHIVIATORE_VOLUME_CLEANER_H

#include <vector>
#include <cstdio>
#include <algorithm>
#include <iostream>
namespace cumbac {

unsigned get_new_cell_num(unsigned orig_cell_num, unsigned max_range=0);

template <class TB>
struct Beams
{
  TB data_z[1024];
  TB data_d[1024];
  TB data_v[1024];
  TB data_w[1024];
};

template <class TB>
struct BeamCleaner
{
  TB bin_wind_magic_number;
  TB Z_missing;
  TB W_threshold;
  TB V_missing;

  unsigned min_segment_length =  4;
  unsigned max_segment_length = 40;

//  BeamCleaner();

//  void print_config(FILE* out) const;

//  void clean_beams(Beams<TB>& b, unsigned beam_size, std::vector<bool>& corrected) const;

	BeamCleaner(TB V, TB Z, TB W)
	{
	  bin_wind_magic_number =  V;
	  Z_missing             =  Z;
	  W_threshold           =  W;
	  V_missing		=  V;
	}
	BeamCleaner(TB V, TB Z, TB W, TB VM)
	{
	  bin_wind_magic_number =  V;
	  Z_missing             =  Z;
	  W_threshold           =  W;
	  V_missing		= VM;
//std::cout<<bin_wind_magic_number<<" "<<Z_missing<<" "<<W_threshold<<" "<<V_missing<<std::endl;
	}


	void print_config(FILE* out) const
	{
	    fprintf(out, "bin_wind_magic_number: %u, min_segment_length: %u, max_segment_length: %u\n",
	      (unsigned)bin_wind_magic_number, min_segment_length, max_segment_length);
	}
	

void clean_beams(Beams<TB>& b, unsigned beam_size, std::vector<bool>& corrected) const
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
unsigned counter = 0;

  for (ibin =0; ibin < beam_size; ibin++)
  {
    if (!in_a_segment) {
/* cerco la prima cella segmento da pulire*/
      //if (b.data_w[ibin] == 0 && b.data_v_w[ibin] == -125 )
  //    std::cout<<ibin<<" "<<b.data_w[ibin]<<" "<< b.data_v[ibin]<<std::endl;
      if (b.data_w[ibin] == W_threshold && b.data_v[ibin] == bin_wind_magic_number ){
        in_a_segment = true;
        start = ibin;
        after = false;
        before = false;
      }
    } else {
/* cerco la fine segmento da pulire*/
      //if (b.data_w[ibin] != 0 || b.data_v_w[ibin] != -125 || ibin == (beam_info.cell_num -1) )
      if (b.data_w[ibin] != W_threshold || b.data_v[ibin] != bin_wind_magic_number || ibin == (beam_size - 1) ){
        in_a_segment = false;
        end = ibin -1;
        if ( ibin == (beam_size -1) ) end=ibin;  // caso particolare per fine raggio
/* Fine trovata ora procedo alla pulizia eventuale */
        segment_length = end - start;
	counter = counter + (unsigned)(segment_length);
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
            b.data_z[ib] = Z_missing;
            b.data_w[ib] = W_threshold;
            b.data_v[ib] = V_missing;
//            *((signed char*)b.data_v + ib) = -128;
            corrected[ib] = true;
          }
        }
      }
    }
  }
//std::cout<<"trovate #"<<counter<<" celle corrette"<<std::endl;
}



};

}

#endif
