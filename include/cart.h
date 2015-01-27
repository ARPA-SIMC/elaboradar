#ifndef ELABORADAR_CART_H
#define ELABORADAR_CART_H

#include <elaboradar/matrix.h>
#include <elaboradar/volume.h>

namespace elaboradar {

class CartFullRes
{
public:
    /// Missing value in the azimuth and range index mappings
    static const unsigned short missing = 0xffff;

    const unsigned beam_size;
    /// Azimuth indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<unsigned short> map_azimuth;
    /// Range indices to use to lookup a map point in a volume
    /// -1 means no mapping
    Matrix2D<unsigned short> map_range;

    /**
     * Build a cartography mapping cartesian coordinates to volume polar
     * indices.
     *
     * The mapping is a 1 to 1 mapping, without scaling.
     *
     * ignore_data is only used for tests: when true, the data in the scan will
     * be ignored, so that tests can test the correctness of the spatial
     * mapping functions regardless of the actual data in the volume.
     */
    CartFullRes(const PolarScan<double>& scan, bool ignore_data=false);

#if 0
    /// Get the value of a polar scan at the point corresponding to these x and
    /// y values in the map
    template<typename T>
    void to_cart(const Matrix2D<T>& src, Matrix2D<T>& dst)
    {
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
#endif
};

}

#endif
