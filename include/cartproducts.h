#ifndef ELABORADAR_CART_PRODUCTS_H
#define ELABORADAR_CART_PRODUCTS_H

#include <radarelab/cart.h>

namespace radarelab {

struct Assets;

struct CartProducts
{
    log4c_category_t* logging_category;

    // Main polar -> cartesian coordinate mapping
    CoordinateMapping mapping;

    // Coordinate mapping for full size images, selecting the point where the
    // volume has max dBZ
    FullsizeIndexMapping fullres;

    // Coordinate mapping for scaled images, selecting the point where the
    // volume has max dBZ
    ScaledIndexMapping scaled;

    Image<unsigned char> z_out;
    Image<unsigned char> top_1x1;
    Image<unsigned char> qual_Z_1x1;
    Image<unsigned char> quota_1x1;
    Image<unsigned char> dato_corr_1x1;
    Image<unsigned char> elev_fin_1x1;
    Image<unsigned char> beam_blocking_1x1;
    Image<unsigned char> neve_1x1;
    Image<unsigned char> corr_1x1;
    Image<unsigned char> conv_1x1;

    /**
     * Hold products for volume, with square images with a side of \a
     * image_side pixels, and for each pixel sampling a square area \a
     * sample_square_size pixels wide
     */
    CartProducts(const Volume<double>& volume, unsigned image_side, unsigned sample_square_size);

    void write_out(Assets& assets);
};

}

#endif
