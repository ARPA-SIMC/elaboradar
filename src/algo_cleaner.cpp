/*
 * =====================================================================================
 *
 *       Filename:  volume_cleaner.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  18/02/2014 12:19:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "algo/cleaner.h"

namespace elaboradar {
namespace algo {

using namespace std;

std::vector<bool> Cleaner::clean_beam(const Eigen::VectorXd& beam_z, const Eigen::VectorXd& beam_w, const Eigen::VectorXd& beam_v) const
{
    const unsigned beam_size = beam_z.rows();
    vector<bool> res(beam_size, false);
    bool in_a_segment = false;
    unsigned start, end;
    unsigned segment_length;
    bool before, after;
    unsigned counter = 0;

    for (unsigned ibin = 0; ibin < beam_size; ++ibin)
    {
        if (!in_a_segment)
        {
            /* cerco la prima cella segmento da pulire*/
            //if (b.data_w[ibin] == 0 && b.data_v_w[ibin] == -125 )
            //    std::cout<<ibin<<" "<<b.data_w[ibin]<<" "<< b.data_v[ibin]<<std::endl;
            if (beam_w(ibin) == W_threshold && beam_v(ibin) == bin_wind_magic_number)
            {
                in_a_segment = true;
                start = ibin;
                after = false;
                before = false;
            }
        } else {
            /* cerco la fine segmento da pulire*/
            //if (b.data_w[ibin] != 0 || b.data_v_w[ibin] != -125 || ibin == (beam_info.cell_num -1) )
            if (beam_w(ibin) != W_threshold || beam_v(ibin) != bin_wind_magic_number || ibin == (beam_size - 1))
            {
                in_a_segment = false;
                end = ibin - 1;
                if (ibin == (beam_size - 1)) end = ibin;  // caso particolare per fine raggio
                /* Fine trovata ora procedo alla pulizia eventuale */
                segment_length = end - start;
                counter = counter + (unsigned)(segment_length);

                /* Cerco dati validi in Z prima del segmento */
                for (int ib = ibin - 12; ib < (signed)ibin; ++ib)
                    if (ib >= 0 && beam_z(ib) > Z_missing)
                        before = true;

                /* Cerco dati validi in Z dopo il segmento */
                for (unsigned ia = ibin + 1; ia <= ibin + 12; ++ia)
                    if (ia < beam_size && beam_z(ia) >= Z_missing)
                        after = true;

                if ((segment_length >= min_segment_length && !before && !after) ||
                        segment_length >= max_segment_length)
                {
                    /* qui pulisco */
                    //         printf (" pulisco %d %d %d \n",segment_length, min_segment_length, max_segment_length);
                    for (unsigned ib = start; ib <= end; ++ib)
                        res[ib] = true;
                }
            }
        }
    }
    return res;
}

void Cleaner::clean(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v)
{
    return clean(scan_z, scan_w, scan_v, scan_v.undetect);
}

void Cleaner::clean(PolarScan<double>& scan_z, PolarScan<double>& scan_w, PolarScan<double>& scan_v, double bin_wind_magic_number)
{
    if (scan_z.beam_count != scan_w.beam_count)
        throw std::runtime_error("scan_z beam_count is different than scan_w beam_count");
    if (scan_z.beam_size != scan_w.beam_size)
        throw std::runtime_error("scan_z beam_size is different than scan_w beam_size");

    if (scan_z.beam_count != scan_v.beam_count)
        throw std::runtime_error("scan_z beam_count is different than scan_v beam_count");
    if (scan_z.beam_size != scan_v.beam_size)
        throw std::runtime_error("scan_z beam_size is different than scan_v beam_size");

    Cleaner cleaner(scan_z.undetect, scan_w.undetect, scan_v.nodata, bin_wind_magic_number);

    const unsigned beam_count = scan_z.beam_count;
    const unsigned beam_size = scan_z.beam_size;

    //fprintf(stderr, "NEWCLEANER zmis %f, wthr %f, vmis %f, mn %f\n",
    //        cleaner.Z_missing, cleaner.W_threshold, cleaner.V_missing, cleaner.bin_wind_magic_number);

    for (unsigned i = 0; i < beam_count; ++i)
    {
        // Compute which elements need to be cleaned
        vector<bool> corrected = cleaner.clean_beam(scan_z.row(i), scan_w.row(i), scan_v.row(i));

        for (unsigned ib = 0; ib < beam_size; ++ib)
            if (corrected[ib])
            {
                scan_z(i, ib) = cleaner.Z_missing;
                scan_w(i, ib) = cleaner.W_threshold;
                scan_v(i, ib) = cleaner.V_missing;
            }
    }
}

}
}
