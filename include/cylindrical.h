#ifndef ARCHIVIATORE_CILINDRICAL_CLASS_H
#define ARCHIVIATORE_CILINDRICAL_CLASS_H

#include <vector>
#include <volume.h>

namespace cumbac {

struct CylindricalVolume
{
    std::vector<Matrix2D<double>*> slices;
    const unsigned x_size;
    const unsigned z_size;

    CylindricalVolume(unsigned slice_count, unsigned x_size, unsigned z_size, double missing_value)
        : x_size(x_size), z_size(z_size)
    {
        slices.reserve(slice_count);
        for (unsigned i = 0; i < slice_count; ++i)
            slices.push_back(new Matrix2D<double>(Matrix2D<double>::Constant(x_size, z_size, missing_value)));
    }
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

    void resample(const Volume<double>& volume, unsigned max_bin, double size_cell);
};

}

#endif
