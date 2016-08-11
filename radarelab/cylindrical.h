/**
 *  @file
 *  @ingroup radarelab 
*/
#ifndef RADARELAB_CILINDRICAL_CLASS_H
#define RADARELAB_CILINDRICAL_CLASS_H

#include <vector>
#include <radarelab/volume.h>

namespace radarelab {

/**
 * Radar volume mapped to cylindrical coordinates.
 */
struct CylindricalVolume
{
    /**
     * Vertical rectangular x,z semi-slices of the cylinder, with one side resting
     * on the cylinder axis. x grows along the radius of the cylinder, z grows
     * along the axis.
     *
     * The angle between slices is constant.
     */
    std::vector<Matrix2D<double>*> slices;
    unsigned x_size;
    unsigned z_size;

    /// Resolution in x and z
    double resol[2];

    CylindricalVolume(const Volume<double>& volume, unsigned slice_count, double missing_value, double x_res, double z_res);
    ~CylindricalVolume()
    {
        for (std::vector<Matrix2D<double>*>::iterator i = slices.begin(); i != slices.end(); ++i)
            delete *i;
    }

    double& operator()(unsigned slice, unsigned row, unsigned col)
    {
        //if (slice >= slices.size()) throw std::runtime_error("slices: fuori coordinata slice");
        return (*slices[slice])(row, col);
    }

    const double& operator()(unsigned slice, unsigned row, unsigned col) const
    {
        //if (slice >= slices.size()) throw std::runtime_error("slices: fuori coordinata slice");
        return (*slices[slice])(row, col);
    }

    void resample(const Volume<double>& volume, unsigned max_bin);
};

}

#endif
