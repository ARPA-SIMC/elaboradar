#include "volume.h"

namespace cumbac {

template<typename T>
void volume_resample(const Volume<T>& src, Volume<T>& dst)
{
}

template void volume_resample<double>(const Volume<double>& src, Volume<double>& dst);
template void volume_resample<unsigned char>(const Volume<unsigned char>& src, Volume<unsigned char>& dst);

}
