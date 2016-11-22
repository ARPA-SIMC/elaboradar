#ifndef ELABORADAR_CART_PRODUCTS_H
#define ELABORADAR_CART_PRODUCTS_H

#include <radarelab/cart.h>

namespace elaboradar {

struct Assets;

/* structure to describe general features for an odim product */ 
struct OdimProdDefs 
{
    radarelab::Image<unsigned char> & prodField;
    std::string System ; 
    std::string ProductType ;
    std::string Quantity ;
    unsigned char Nodata ;
    unsigned char Undetect;
    double Offset ;
    double QuantityDynamics ;
    bool SaveQuality;
    radarelab::Image<unsigned char> & QualityField;
    double prodRes;
    std::vector<double> prodPar;

	OdimProdDefs(radarelab::Image<unsigned char> & prodField, double prodRes) ;
	OdimProdDefs(radarelab::Image<unsigned char> & prodField, radarelab::Image<unsigned char> & QualityField, double prodRes) ;

};

struct CartProducts
{
    log4c_category_t* logging_category;

    // Main polar -> cartesian coordinate mapping
    radarelab::CoordinateMapping mapping;

    // Coordinate mapping for full size images, selecting the point where the
    // volume has max dBZ
    radarelab::FullsizeIndexMapping fullres;

    // Coordinate mapping for scaled images, selecting the point where the
    // volume has max dBZ
    radarelab::ScaledIndexMapping scaled;

    radarelab::Image<unsigned char> z_out;
    radarelab::Image<unsigned char> z_fr;
    radarelab::Image<unsigned char> top_1x1;
    radarelab::Image<unsigned char> qual_Z_1x1;
    radarelab::Image<unsigned char> quota_1x1;
    radarelab::Image<unsigned char> dato_corr_1x1;
    radarelab::Image<unsigned char> elev_fin_1x1;
    radarelab::Image<unsigned char> beam_blocking_1x1;
    radarelab::Image<unsigned char> neve_1x1;
    radarelab::Image<unsigned char> corr_1x1;
    radarelab::Image<unsigned char> conv_1x1;

	double FullsizeRes = 0.;
	double ScaledRes   = 0.;

    /**
     * Hold products for volume, with square images with a side of \a
     * image_side pixels, and for each pixel sampling a square area \a
     * sample_square_size pixels wide
     */
    CartProducts(const radarelab::Volume<double>& volume, unsigned image_side, unsigned sample_square_size);

    void write_out(Assets& assets);
    void write_out(Assets& assets, unsigned image_side);
    void write_out(Assets& assets, unsigned image_side, std::string algos);

    void write_odim(Assets& assets, unsigned image_side, std::string algos, OdimProdDefs & odimProd);

};


}

#endif
