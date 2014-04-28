#ifndef ARCHIVIATORE_VOLUME_ELEV_FIN_CLASS_H
#define ARCHIVIATORE_VOLUME_ELEV_FIN_CLASS_H

#include "volume.h"
#include "volume/loader.h"
#include <H5Cpp.h>

namespace cumbac {
namespace volume {

template<typename T>
struct ElevFin
{
    const Volume<T>& volume;
    const volume::LoadInfo& load_info;

    // elevazione finale in coordinate azimut range
    std::vector<unsigned char> elev_fin[NUM_AZ_X_PPI];

    ElevFin(const Volume<T>& volume, const volume::LoadInfo& load_info)
        : volume(volume), load_info(load_info) {}

    std::vector<unsigned char>& operator[](unsigned idx) { return elev_fin[idx]; }
    const std::vector<unsigned char>& operator[](unsigned idx) const { return elev_fin[idx]; }

    void init()
    {
        // FIXME: set to 0 to have the right size. We start from 512 (MAX_BIN)
        // to allocate enough memory for legacy code that iterates on MAX_BIN
        // to successfully read zeroes
        unsigned max_size = 512;
        for (unsigned iel = 0; iel < volume.size(); ++iel)
        {
            if (volume.scan(iel).beam_size && volume.scan(iel).beam_size > max_size)
                max_size = volume.scan(iel).beam_size;
        }

        for (unsigned i = 0; i < NUM_AZ_X_PPI; ++i)
        {
            elev_fin[i].resize(max_size, 0);
        }
    }

    inline double elevation_rad_at_elev_preci(unsigned az_idx, unsigned ray_idx) const
    {
        unsigned az = elev_fin[az_idx][ray_idx];
        return load_info.scan(az).get_elevation_rad(az_idx);
    }

    inline double elevation_at_elev_preci(unsigned az_idx, unsigned ray_idx) const
    {
        unsigned az = elev_fin[az_idx][ray_idx];
        return load_info.scan(az).get_elevation(az_idx);
    }

    inline double db_at_elev_preci(unsigned az_idx, unsigned ray_idx) const
    {
        const PolarScan<T>& s = volume.scan(elev_fin[az_idx][ray_idx]);
        if (ray_idx < s.beam_size)
            return s.get(az_idx, ray_idx);
        else
            // If we are reading out of bounds, return 1 (the missing value)
            return MINVAL_DB;
    }

    void write_info_to_debug_file(H5::H5File out)
    {
        using namespace H5;

        // Compute dimensions
        hsize_t dims[2] = { NUM_AZ_X_PPI, 0 };
        for (unsigned i = 0; i < NUM_AZ_X_PPI; ++i)
            if (dims[1] < elev_fin[i].size())
                dims[1] = elev_fin[i].size();

        DataSpace file_data_space(2, dims);

        // Dataset data type
        IntType datatype( PredType::NATIVE_UCHAR );

        // Dataset fill value
        DSetCreatPropList props;
        unsigned char fill_value(0);
        props.setFillValue(datatype, &fill_value);

        // Create the dataset
        DataSet ds = out.createDataSet("/elev_fin", datatype, file_data_space, props);

        // Write elev_fin to it
        for (unsigned i = 0; i < NUM_AZ_X_PPI; ++i)
        {
            hsize_t mdims[1] = { elev_fin[i].size() };
            DataSpace memory_data_space(1, mdims);

            hsize_t count[] = { 1, elev_fin[i].size() };
            hsize_t start[] = { i, 0 };
            file_data_space.selectHyperslab(H5S_SELECT_SET, count, start);

            ds.write(elev_fin[i].data(), datatype, memory_data_space, file_data_space);
        }
    }
};

}
}

#endif
