#ifndef ELABORADAR_CART_H
#define ELABORADAR_CART_H

#include <elaboradar/matrix.h>
#include <elaboradar/volume.h>
#include <limits>

namespace elaboradar {

/**
 * Mapping of cartesian coordinates to raw azimuth angles and range distances.
 */
struct CoordinateMapping
{
    const unsigned beam_size;
    /// Azimuth indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<double> map_azimuth;
    /// Range indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<double> map_range;

    /**
     * Build a cartography mapping cartesian coordinates to volume polar
     * indices.
     *
     * The mapping is a 1 to 1 mapping, without scaling.
     */
    CoordinateMapping(unsigned beam_size);
};


/**
 * Mapping of cartesian coordinates to specific azimuth and range volume indices
 */
class IndexMapping
{
public:
    /// Missing value in the azimuth and range index mappings
    static const unsigned missing = 0xffffffff;

    const unsigned beam_size;
    /// Azimuth indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<unsigned> map_azimuth;
    /// Range indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<unsigned> map_range;

    IndexMapping(unsigned beam_size);

    /**
     * Map cartesian cardinates to polar volume indices. When a cartesian
     * coordinate maps to more than one polar value, take the one with the
     * maximum data value.
     */
    void map_max_sample(const PolarScan<double>& scan);

    /**
     * Same as map_max_sample(PolarScan), but reuse an existing
     * CoordinateMapping.
     */
    void map_max_sample(const PolarScan<double>& scan, const CoordinateMapping& mapping);

    /// Copy data from the polar scan src to the cartesian map dst
    template<typename T>
    void to_cart(const PolarScan<T>& src, Matrix2D<T>& dst)
    {
        // In case dst is not a square with side beam_size*2, center it
        int dx = (beam_size * 2 + dst.cols()) / 2;
        int dy = (beam_size * 2 + dst.rows()) / 2;

        for (unsigned y = 0; y < dst.rows(); ++y)
        {
            if (y + dy < 0 || y + dy >= beam_size * 2) continue;

            for (unsigned x = 0; x < dst.cols(); ++x)
            {
                if (x + dx < 0 || x + dx >= beam_size * 2) continue;

                auto azimuth = map_azimuth(y + dy, x + dx);
                auto range = map_range(y + dy, x + dx);

                if (azimuth == missing || range == missing) continue;
                if (azimuth >= src.beam_count || range >= src.beam_size) continue;

                dst(y, x) = src(azimuth, range);
            }
        }
    }
};

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

#endif
