#include "cart.h"
#include <cmath>

namespace elaboradar {

CoordinateMapping::CoordinateMapping(unsigned beam_size)
    : beam_size(beam_size),
      map_azimuth(beam_size * 2, beam_size * 2),
      map_range(beam_size * 2, beam_size * 2)
{
    for (unsigned y = 0; y < beam_size * 2; ++y)
        for (unsigned x = 0; x < beam_size * 2; ++x)
        {
            // x and y centered on the map center
            double absx = x + 0.5 - beam_size;
            double absy = y + 0.5 - beam_size;

            // Compute range
            map_range(y, x) = hypot(absx, absy);

            // Compute azimuth
            double az;
            if (absx > 0)
                if (absy < 0)
                    az = atan(absx / -absy) * M_1_PI * 180.;
                else
                    az = 90.0 + atan(absy / absx) * M_1_PI * 180.;
            else
                if (absy > 0)
                    az = 180 + atan(-absx / absy) * M_1_PI * 180.;
                else
                    az = 270 + atan(-absy / -absx) * M_1_PI * 180.;
            map_azimuth(y, x) = az;
        }
}


const unsigned IndexMapping::missing;

IndexMapping::IndexMapping(unsigned height, unsigned width)
    : height(height), width(width),
      map_azimuth(Matrix2D<unsigned>::Constant(height, width, missing)),
      map_range(Matrix2D<unsigned>::Constant(height, width, missing))
{
}

FullsizeIndexMapping::FullsizeIndexMapping(unsigned beam_size)
    : IndexMapping(beam_size * 2, beam_size * 2)
{
}

void FullsizeIndexMapping::map_max_sample(const PolarScan<double>& scan)
{
    CoordinateMapping raw_mapping(scan.beam_size);
    map_max_sample(scan, raw_mapping);
}

void FullsizeIndexMapping::map_max_sample(const PolarScan<double>& scan, const CoordinateMapping& mapping)
{
    for (unsigned y = 0; y < height; ++y)
        for (unsigned x = 0; x < width; ++x)
        {
            // Map cartesian coordinates to angles and distances
            unsigned range = floor(mapping.map_range(y, x));
            if (range >= scan.beam_size) continue;

            // Exact angle
            double az = mapping.map_azimuth(y, x);

            // Iterate indices 0.45° before and after
            int az_min = floor((az - .45) * scan.beam_count / 360);
            int az_max = ceil((az + .45) * scan.beam_count / 360);
            if (az_min < 0)
            {
                az_min += scan.beam_count;
                az_max += scan.beam_count;
            }

            // Look for the position of the maximum value of the scan in
            // this range
            bool first = true;
            double maxval;
            unsigned maxval_idx;
            for (unsigned iaz = az_min; iaz < (unsigned)az_max; ++iaz)
            {
                unsigned azidx = iaz % scan.beam_count;
                unsigned char sample = scan.get(azidx, range);

                if (first || sample > maxval)
                {
                    maxval = sample;
                    maxval_idx = azidx;
                    first = false;
                }
            }

            // If nothing has been found, skip this sample
            if (first) continue;

            // Store the indices for this mapping
            map_azimuth(y, x) = maxval_idx;
            map_range(y, x) = range;
        }
}


ScaledIndexMapping::ScaledIndexMapping(unsigned beam_size, unsigned image_side, unsigned fullsize_pixels_per_scaled_pixel)
    : IndexMapping(image_side, image_side),
      fullsize_pixels_per_scaled_pixel(fullsize_pixels_per_scaled_pixel),
      image_offset((beam_size * 2 - image_side * fullsize_pixels_per_scaled_pixel) / 2)
{
    if ((image_side * fullsize_pixels_per_scaled_pixel) % 2 != 0)
        throw std::runtime_error("the image cannot be properly centered on the full size image");
}

void ScaledIndexMapping::map_max_sample(const PolarScan<double>& scan, const FullsizeIndexMapping& mapping)
{
    // ciclo sui punti della nuova matrice. per il primo prenderò il massimo tra i primi sedici etc..
    for (unsigned y = 0; y < height; ++y)
        for (unsigned x = 0; x < width; ++x)
        {
            bool first = true;
            double maxval;
            unsigned maxval_az = missing;
            unsigned maxval_range = missing;

            // Load each sample with a value from a 4x4 window on the original image
            for(unsigned sy = 0; sy < fullsize_pixels_per_scaled_pixel; ++sy)
                for(unsigned sx = 0; sx < fullsize_pixels_per_scaled_pixel; ++sx)
                {
                    // Use the full size mapping to get the volume value at this point
                    unsigned src_x = x * fullsize_pixels_per_scaled_pixel + sx + image_offset;
                    unsigned src_y = y * fullsize_pixels_per_scaled_pixel + sy + image_offset;
                    unsigned range = mapping.map_range(src_y, src_x);
                    if (range == missing) continue;
                    unsigned az = mapping.map_azimuth(src_y, src_x);
                    unsigned char sample = scan.get(az, range);

                    // Find the source point with the maximum sample
                    if (first || sample > maxval)
                    {
                        maxval = sample;
                        maxval_az = az;
                        maxval_range = range;
                        first = false;
                    }
                }

            // If nothing has been found, skip this sample
            if (first) continue;

            // Store the indices for this mapping
            map_azimuth(y, x) = maxval_az;
            map_range(y, x) = maxval_range;
        }
}

#if 0
/**
 * Scaled version of CartFullres
 */
class CartScaled
{
    const unsigned beam_size;
    /// Azimuth indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<unsigned> map_azimuth;
    /// Range indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<unsigned> map_range;

        for(unsigned quad = 0; quad < 4; ++quad)
            for(unsigned qy = 0; qy < beam_size; ++qy)
                for(unsigned qx = 0; qx < beam_size; ++qx)
                {
                    unsigned x;
                    unsigned y;
                    double az;
                    int irange = map_range(qy, qx);
                    if (irange == -1) continue;
                    switch(quad)
                    {
                        case 0:
                            x = beam_size + qx;
                            y = beam_size - qy;
                            az = map_azimut(qy, qx);
                            break;
                        case 1:
                            x = beam_size + qx;
                            y = beam_size + qy;
                            az = map_azimut(qy, qx) + 90.;
                            break;
                        case 2:
                            x = beam_size - qx;
                            y = beam_size + qy;
                            az = map_azimut(qy, qx) + 180.;
                            break;
                        case 3:
                            x = beam_size - qx;
                            y = beam_size - qy;
                            az = map_azimut(qy, qx) + 270.;
                            break;
                    }

                    int az_min = (int)((az - .45)/.9);
                    int az_max = ceil((az + .45)/.9);


                    if (az_min < 0)
                    {
                        az_min = az_min + NUM_AZ_X_PPI;
                        az_max = az_max + NUM_AZ_X_PPI;
                    }

                    cont=0;
                    for (unsigned iaz = az_min; iaz < az_max; ++iaz)
                    {
                        // Enrico: cerca di non leggere fuori dal volume effettivo
                        unsigned char sample = 0;
                        if (irange < cb.volume[0].beam_size)
                            sample = max(DBtoBYTE(cb.volume[0].get(iaz%NUM_AZ_X_PPI, irange)), (unsigned char)1);   // il max serve perchè il valore di MISSING è 0

                        if(cart(y, x) <= sample){
                            cart(y, x) = sample;
                            topxy(y, x)=cb.top(iaz%NUM_AZ_X_PPI, irange);
                            if (cb.do_quality)
                            {
                                if (irange < cb.volume[cb.anaprop.elev_fin[iaz%NUM_AZ_X_PPI][irange]].beam_size)
                                    qual_Z_cart(y, x) = cb.qual->scan(cb.anaprop.elev_fin[iaz%NUM_AZ_X_PPI][irange]).get(iaz%NUM_AZ_X_PPI, irange);
                                else
                                    qual_Z_cart(y, x) = 0;
                                quota_cart(y, x)=cb.anaprop.quota(iaz%NUM_AZ_X_PPI, irange);
                                // if (iaz == 0) LOG_DEBUG(" x,y %4d,%4d - irange %4d quota %d",x,y,irange,cb.anaprop.quota(iaz%NUM_AZ_X_PPI, irange));
                                dato_corr_xy(y, x)=cb.anaprop.dato_corrotto(iaz%NUM_AZ_X_PPI, irange);
                                beam_blocking_xy(y, x)=cb.beam_blocking(iaz%NUM_AZ_X_PPI, irange);
                                elev_fin_xy(y, x)=cb.anaprop.elev_fin[iaz%NUM_AZ_X_PPI][irange];
                                /*neve_cart(y, x)=qual_Z_cart(y, x);*/
                                if (cb.calcolo_vpr)
                                {
                                    neve_cart(y, x)=(cb.calcolo_vpr->neve(iaz%NUM_AZ_X_PPI, irange))?0:1;
                                    corr_cart(y, x)=cb.calcolo_vpr->corr_polar(iaz%NUM_AZ_X_PPI, irange);
                                }
                            }
                            if (cb.do_class)
                            {
                                if (irange<cb.calcolo_vpr->x_size)
                                    conv_cart(y, x)=cb.calcolo_vpr->conv(iaz%NUM_AZ_X_PPI,irange);
                            }
                        }
                        if (cb.do_zlr_media)
                        {
                            //if (volume.scan(0).get_raw(iaz%NUM_AZ_X_PPI, irange) > 0)
                            if (sample > 0)
                            {
                                cartm(y, x)=cartm(y, x)+BYTEtoZ(sample);
                                cont=cont+1;
                            }
                        }
                    }

                    if (cb.do_zlr_media)
                    {
                        if (cont > 0) cartm(y, x)=cartm(y, x)/(float)(cont);
                    }
                    /*
                     *****  per scrivere in griglia cartesiana************
                     bloc_xy[MAX_BIN*2-y][x]=beam_blocking((int)((float)(az)/.9), irange);
                     elev_xy[MAX_BIN*2-y][x]=volume.elev_fin[(int)((float)(az)/.9)][irange];
                     dato_corrotto_xy(MAX_BIN*2-y, x)= dato_corrotto((int)((float)(az)/.9), irange);
                     elev_fin_xy[MAX_BIN*2-y][x]=first_level((int)((float)(az)/.9), irange);
                     dato_corrotto_xy(MAX_BIN*2-y, x)= dato_corrotto((int)((float)(az)/.9), irange); */
                }
    }
};
#endif

}
