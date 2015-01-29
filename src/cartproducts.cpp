#include "cartproducts.h"
#include "assets.h"

namespace elaboradar {

CartProducts::CartProducts(const Volume<double>& volume, unsigned image_side, unsigned sample_square_size)
    : fullres(volume[0].beam_size), scaled(volume[0].beam_size, image_side, sample_square_size),
      z_out(image_side), top_1x1(image_side), qual_Z_1x1(image_side),
      quota_1x1(image_side), dato_corr_1x1(image_side),
      elev_fin_1x1(image_side), beam_blocking_1x1(image_side),
      neve_1x1(image_side), corr_1x1(image_side), conv_1x1(image_side)
{
    logging_category = log4c_category_get("radar.cart");

    quota_1x1.fill(128);

    LOG_INFO("Creazione Matrice Cartesiana");
    fullres.map_max_sample(volume[0]);
    //assets.write_gdal_image(fullres.map_azimuth, "DIR_DEBUG", "map_azimuth", "PNG");
    //assets.write_gdal_image(fullres.map_range, "DIR_DEBUG", "map_range", "PNG");

    LOG_INFO("Creazione Matrice Cartesiana ridimensionata");
    scaled.map_max_sample(volume[0], fullres);
}

void CartProducts::write_out(Assets& assets)
{
    assets.write_image(z_out, "OUTPUT_Z_LOWRIS_DIR", ".ZLR", "file output 1X1");
    assets.write_image(top_1x1, "DIR_QUALITY", ".top20_ZLR", "file top20");
    assets.write_image(qual_Z_1x1, "OUTPUT_Z_LOWRIS_DIR", ".qual_ZLR", "file qualita' Z");
    assets.write_image(quota_1x1, "DIR_QUALITY", ".quota_ZLR", "file qel1uota");
    assets.write_image(dato_corr_1x1, "DIR_QUALITY", ".anap_ZLR", "file anap");
    assets.write_image(elev_fin_1x1, "DIR_QUALITY", ".elev_ZLR", "file elev");
    assets.write_image(beam_blocking_1x1, "DIR_QUALITY", ".bloc_ZLR", "file bloc");
    assets.write_image(neve_1x1, "DIR_QUALITY", ".corr_ZLR", "file correzione VPR");
    assets.write_image(corr_1x1, "DIR_QUALITY", ".corr_ZLR", "file correzione VPR");
    assets.write_image(conv_1x1, "DIR_QUALITY", ".conv_ZLR", "punti convettivi");
}

}
