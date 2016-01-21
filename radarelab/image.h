/**
 *  @file
 *  @ingroup radarelab 
*/
#ifndef RADARELAB_IMAGE_H
#define RADARELAB_IMAGE_H

#include <radarelab/matrix.h>
#include <stdexcept>
#include <string>

namespace radarelab {

/**
 * Initialize the GDAL library when called for the first time; does nothing all
 * other times.
 */
void gdal_init_once();

template<typename T>
void write_image(const Matrix2D<T>& image, const std::string& fname, const std::string& format);

std::string gdal_extension_for_format(const std::string& format);

}

#endif
