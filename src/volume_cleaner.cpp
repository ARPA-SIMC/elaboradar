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
#include "volume_cleaner.h"

namespace cumbac {

using namespace std;

// TODO: toglierlo
#define MAX_DIM 512

// TODO: direi che nessuno usa questa funzione, si potrebbe cassare
unsigned get_new_cell_num(unsigned orig_cell_num, unsigned max_range)
{
    // Lunghezza che vogliamo
    unsigned rmax = min(orig_cell_num, (unsigned)MAX_DIM);
    if (max_range > 0)
      rmax = min(rmax, max_range);
    return rmax;
}


//template struct BeamCleaner <unsigned char>;  // perchè ci va messo?
template struct BeamCleaner <double>;  // perchè ci va messo?


}
