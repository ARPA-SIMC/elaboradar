#ifndef ELABORADAR_ALGO_UTILS_H
#define ELABORADAR_ALGO_UTILS_H

#include <math.h>

namespace elaboradar {
namespace algo {

/**
 * Compute the corrected value (in dB) given the original value (in dB) and the
 * beam blocking percentage (from 0 to 100) for that value
 */
static inline double beam_blocking_correction(double val_db, double beamblocking)
{
   return val_db - 10 * log10(1. - beamblocking / 100.);
}



}
}

#endif
