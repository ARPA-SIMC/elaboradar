#include "volume.h"
#include "logging.h"
#include <stdexcept>
#include <algorithm>

//#define IMPRECISE_AZIMUT

using namespace std;

namespace cumbac {

void VolumeStats::print(FILE* out)
{
    fprintf(out, "Nel    Zeros     Ones   Others      Sum\n");
    for (size_t iel =0; iel<count_zeros.size(); ++iel){
        fprintf(out, "%4zu %8u %8u %8u %8u\n",iel,count_zeros[iel],count_ones[iel],count_others[iel],sum_others[iel]);
    }
}

template class Volume<double>;

}
