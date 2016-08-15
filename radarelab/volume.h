/**
 *  @file
 *  @ingroup radarelab
 *  @brief Definisce le principali strutture che contengono i dati
 */
#ifndef RADARELAB_VOLUME_H
#define RADARELAB_VOLUME_H

#include <radarelab/logging.h>
#include <radarelab/matrix.h>
#include <radarelab/algo/dbz.h>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <memory>
#include <Eigen/Core>
#include <radarelab/RadarSite.h>

// TODO: prima o poi arriviamo a far senza di questi define
#define NUM_AZ_X_PPI 400

#define MISSING_DB (-20.)
#define MINVAL_DB (80./255.-20.)
#define MAXVAL_DB (60.)

namespace radarelab {

/**
 * Basic structure to describe a polar scan, independently of the type of its
 * samples.
 */
struct PolarScanBase
{
    /// Count of beams in this scan
    unsigned beam_count = 0;

    /// Number of samples in each beam
    unsigned beam_size = 0;

    /// Vector of actual azimuths for each beam
    Eigen::VectorXd azimuths_real;

    /**
     * Nominal elevation of this PolarScan, which may be different from the
     * effective elevation of each single beam
     */
    double elevation = 0;

    /// Vector of actual elevations for each beam
    Eigen::VectorXd elevations_real;

    /// Size of a beam cell in meters
    double cell_size = 0;

    /*
     * Constructor
     * Create a PolarScanBase defining beam_count and beam_size
     * @param [in] beam_count
     * @param [in] beam_size
     */
    PolarScanBase(unsigned beam_count, unsigned beam_size);

    /*
     * Constructor
     * Create a copy of a PolarScanBase
     * @param [in] s  - PolarScanBase to be copied
     */
    PolarScanBase(const PolarScanBase& s);

    /**
     * Height in kilometers (legacy) at range gate for beam elevation  + beam_half_width
     * @param [in] rg - index of range gate  
     * @param [in] beam_half_width - semi-amplitude of the beam
     */
    double height(unsigned rg, double beam_half_width=0.0);

    /**
     *  Height difference in kilometers (legacy) between two range gates
     *  @param [in] rg_start - first range gate index
     *  @param [in] rg_end - second range gate index
     * @return Height difference [km]
     */
    double diff_height(unsigned rg_start, unsigned rg_end);

    /**
     * Return the height (in meters) of the sample at the given cell indexa
     * @param [in] cell_idx - index of range gatea
     * @return - height of cell_idx at PolarScan elevation
     */
    double sample_height(unsigned cell_idx) const;

    /**
     * Return the height of a sample (in meters) given center beam elevation
     * (in degrees), range (in meters) and equivalent earth radius (in meters)
     * @param [in] elevation - Elevation in degrees
     * @param [in] range - range of range gate [m]
     * @param [in] equiv_earth_radius - equivalent earth radius (in meters)
     * @return - height of cell_idx at PolarScan elevation
     */
    static double sample_height(double elevation, double range, double equiv_earth_radius);

    /**
     * Return the height of a sample (in meters) given center beam elevation
     * (in degrees) and range (in meters), using the standard 4/3 equivalent
     * earth radius (in meters)
     * @param [in] elevation - Elevation in degrees
     * @param [in] range - range of range gate [m]
     * @return - height of cell_idx at PolarScan elevation
     */
    static double sample_height(double elevation, double range);
};

/**
 * PolarScan - structure to describe a polarScan containing a matrix of data and conversion factors
 */
template<typename T>
class PolarScan : public PolarScanBase, public Matrix2D<T>
{
public:
    /// Value used as 'no data' value
    T nodata = 0;
    /// Minimum amount that can be measured
    T undetect = 0;
    /// Conversion factor
    T gain = 1;
    /// Conversion factor
    T offset = 0;

    /*!
     * Constructor - Create a polarscan given beam_count, beam_size and the value to fill the Matrix
     * @param [in] beam_count
     * @param [in] beam_size
     * @param[in] default_value
     */
    PolarScan(unsigned beam_count, unsigned beam_size, const T& default_value=algo::DBZ::BYTEtoDB(1))
        : PolarScanBase(beam_count, beam_size),
          Matrix2D<T>(PolarScan::Constant(beam_count, beam_size, default_value))
    {
    }

    /**
     * Constructor
     * Create a copy of a PolarScan
     * @param [in] s  - PolarScan to be copied
     */
    PolarScan(const PolarScan& s) = default;

    template<class OT>
    PolarScan(const PolarScan<OT>& s, const T& default_value)
        : PolarScanBase(s), Matrix2D<T>(PolarScan::Constant(s.beam_count, s.beam_size, default_value)),
          nodata(s.nodata),undetect(s.undetect),gain(s.gain),offset(s.offset)
    {
    }

    ~PolarScan()
    {
    }

    /// Get a beam value
    /// @param az - azimuth index [0-beam_count -1]
    /// @param beam - cell index
    /// @return - bin value
    T get(unsigned az, unsigned beam) const
    {
        return (*this)(az, beam);
    }

    /// Set a beam value
    /// @param az - azimuth index [0-beam_count -1]
    /// @param beam - cell index
    /// @param [in] val - value to be set
    void set(unsigned az, unsigned beam, T val)
    {
        (*this)(az, beam) = val;
    }

    /**
     * Return the height (in meters) of a sample given its azimuth and cell indices
     * use the real beam elevation (not the nominal of the PolarScan) and add half cell_size to the range
     * @param [in] az - azimuth index
     * @param [in] cell_idx - cell index
     * @return height [m]
     */
    double sample_height_real(unsigned az, unsigned cell_idx) const
    {
        return PolarScanBase::sample_height(elevations_real(az), (cell_idx + 0.5) * cell_size);
    }

    /**
     * Fill an array with beam data . If the array is longer than the beam fill the remaining with missing 
     *  @param [in] az - azimuth index
     *  @param [in,out] out - array to be filled
     *  @param out_size - dimension of the array
     *  @param [in] missing - Value to be used to fill the exceeding part 
     */
    void read_beam(unsigned az, T* out, unsigned out_size, T missing=0) const
    {
        using namespace std;

        // Prima riempio il minimo tra ray.size() e out_size
        size_t set_count = min(beam_size, out_size);

        for (unsigned i = 0; i < set_count; ++i)
            out[i] = get(az, i);

        for (unsigned i = set_count; i < out_size; ++i)
            out[i] = missing;
    }

    /**
     * Enlarges the PolarScan increasing beam_size and propagating the last bin value
     * @param [in] new_beam_size 
     */
    void resize_beams_and_propagate_last_bin(unsigned new_beam_size)
    {
        if (new_beam_size <= beam_size) return;
        this->conservativeResize(Eigen::NoChange, new_beam_size);
        this->rightCols(new_beam_size - this->beam_size).colwise() = this->col(this->beam_size - 1);
        this->beam_size = new_beam_size;
    }
};


struct VolumeStats
{
    std::vector<unsigned> count_zeros;
    std::vector<unsigned> count_ones;
    std::vector<unsigned> count_others;
    std::vector<unsigned> sum_others;

    void print(FILE* out);
};

/**
 * Namespace per volume dati
 */
namespace volume {

/**
 * LoadInfo structure - Contains generic volume information 
 */
struct LoadInfo
{
    /// Original file name
    std::string filename;
    // std::string date;
    // std::string time;
    /// Acquisition date
    time_t acq_date;
    /// flag true if data have been decluttered with Doppler at rsp level
    bool declutter_rsp; 

    LoadInfo()
        : declutter_rsp(false)
    {
    }
};

/**
 * Sequence of PolarScans which can have a different beam count for each elevation
 */
template<typename T>
class Scans : public std::vector<PolarScan<T>>
{
public:
    typedef T Scalar;
    /// Odim quantity name
    std::string quantity;
    /// Data units according to ODIM documentation
    std::string units;
    /// Polar volume information
    std::shared_ptr<LoadInfo> load_info;
    /// RadarSite
    RadarSite radarSite;

    Scans() = default;

    /**
     *  Constructor
     *  Copy from another Scans
     *  @param [in] v - Scans object to be copied
     *  @param[in] default_value - default to be used 
     */
    template<typename OT>
    Scans(const Scans<OT>& v, const T& default_value)
    {
        this->quantity = v.quantity;
        this->units = v.units;
        this->load_info = v.load_info;
        this->reserve(v.size());
        this->radarSite = v.radarSite;
        for (const auto& src_scan : v)
            this->push_back(PolarScan<T>(src_scan, default_value));
    }

    /// Access a polar scan
    /// @param [in] idx - index of the PolarScan to be used
    /// @return a pointer to the PolarScan
    PolarScan<T>& scan(unsigned idx) { return (*this)[idx]; }
    /// Access a polar scan (const)
    /// @param [in] idx - index of the PolarScan to be used
    /// @return a pointer to the PolarScan (const)
    const PolarScan<T>& scan(unsigned idx) const { return (*this)[idx]; }

    /**
     * Append a scan to this volume.
     *
     * It is required that scans are added in increasing elevation order,
     * because higher scan indices need to correspond to higher elevation
     * angles.
     *
     * It is required that beam_size is lower than 
     * @param [in] beam_count 
     * @param [in] beam_size
     * @param [in] elevation - PolarScan elevation (degrees)
     * @param [in] cell_size - PolarScan cell size [m]
     */
    PolarScan<T>& append_scan(unsigned beam_count, unsigned beam_size, double elevation, double cell_size)
    {
        // Ensure elevations grow as scan indices grow
        if (!this->empty() && elevation <= this->back().elevation)
        {
            LOG_CATEGORY("radar.io");
            LOG_ERROR("append_scan(beam_count=%u, beam_size=%u, elevation=%f, cell_size=%f) called with an elevation that is not above the last one (%f)", beam_count, beam_size, elevation, cell_size, this->back().elevation);
            throw std::runtime_error("elevation not greather than the last one");
        }

        // Add the new polar scan
        this->push_back(PolarScan<T>(beam_count, beam_size));
        this->back().elevation = elevation;
        this->back().cell_size = cell_size;
        return this->back();
    }

    /**
     * Create or reuse a scan at position idx, with the given beam size 
     * @param [in] idx - index of the PolarScan 
     * @param [in] beam_count 
     * @param [in] beam_size
     * @param [in] elevation
     * @param [in] cell_size
     * @return the idx PolarScan
     */
    PolarScan<T>& make_scan(unsigned idx, unsigned beam_count, unsigned beam_size, double elevation, double cell_size)
    {
        if (idx < this->size())
        {
            if (beam_count != (*this)[idx].beam_count)
            {
                LOG_CATEGORY("radar.io");
                LOG_ERROR("make_scan(idx=%u, beam_count=%u, beam_size=%u) called, but the scan already existed with beam_count=%u", idx, beam_count, beam_size, (*this)[idx].beam_count);
                throw std::runtime_error("beam_size mismatch");
            }
            if (beam_size != (*this)[idx].beam_size)
            {
                LOG_CATEGORY("radar.io");
                LOG_ERROR("make_scan(idx=%u, beam_count=%u, beam_size=%u) called, but the scan already existed with beam_size=%u", idx, beam_count, beam_size, (*this)[idx].beam_size);
                throw std::runtime_error("beam_size mismatch");
            }
        } else {
            // If some elevation has been skipped, fill in the gap
            if (idx > this->size())
            {
                if (this->empty())
                    this->push_back(PolarScan<T>(beam_count, beam_size));
                while (this->size() < idx)
                    this->push_back(PolarScan<T>(beam_count, this->back().beam_size));
            }

            // Add the new polar scan
            this->push_back(PolarScan<T>(beam_count, beam_size));
            this->back().elevation = elevation;
            this->back().cell_size = cell_size;
        }

        // Return it
        return (*this)[idx];
    }

    /**
     * Change the elevations in the PolarScans to match the given elevation vector
     * @param [in] elevations - values to be used
     */
    void normalize_elevations(const std::vector<double>& elevations)
    {
        // Ensure that we have enough standard elevations
        if (elevations.size() < this->size())
        {
            LOG_CATEGORY("radar.io");
            LOG_ERROR("normalize_elevations: standard elevation array has %zd elements, but we have %zd scans", 
                    elevations.size(), this->size());
            throw std::runtime_error("not enough standard elevations");
        }
        // Ensure that the nudging that we do do not confuse a scan
        // with another
        for (size_t i = 0; i < this->size() - 1; ++i)
        {
            if (abs(elevations[i] - this->at(i).elevation) > abs(elevations[i] - this->at(i + 1).elevation))
            {
                LOG_CATEGORY("radar.io");
                LOG_ERROR("normalize_elevations: elevation %zd (%f) should be set to %f but it would make it closer to the next elevation %f", i, this->at(i).elevation, elevations[i], this->at(i + 1).elevation);
                throw std::runtime_error("real elevation is too different than standard elevations");
            }
        }
        // Assign the new elevations
        for (size_t i = 0; i < this->size(); ++i)
            this->at(i).elevation = elevations[i];
    }
};

}

/**
 * Homogeneous volume with a common beam count for all PolarScans
 */
template<typename T>
class Volume : public volume::Scans<T>
{
public:
    /// Number of beam_count used ast each elevations
    const unsigned beam_count;

    /**
     * Constructor 
     * @param [in] beam_count - number of beam_count to be used
     */
    Volume(unsigned beam_count=NUM_AZ_X_PPI)
        : beam_count(beam_count)
    {
    }

    /**
     * Copy constructor
     * @param [in] v - Volume to be copied
     * @param [in] default_value
     */
    template<typename OT>
    Volume(const Volume<OT>& v, const T& default_value)
        : volume::Scans<T>(v, default_value), beam_count(v.beam_count)
    {
    }

    /// Return the maximum beam size in all PolarScans
    const unsigned max_beam_size() const
    {
        unsigned res = 0;
        for (const auto& scan : *this)
            res = std::max(res, scan.beam_size);
        return res;
    }

    /// Test if same cell_size in all PolarScans
    bool is_unique_cell_size() const
    {
        double cell_size = this->at(0).cell_size;
        for (size_t i = 1; i < this->size(); ++i)
            if ( this->at(i).cell_size != cell_size) return false;
        return true;
    }

    /// Return the lowest elevation
    double elevation_min() const
    {
        return this->front().elevation;
    }

    /// Return the highest elevation
    double elevation_max() const
    {
        return this->back().elevation;
    }

    /**
     * Fill a matrix elevations x beam_size with the vertical slice at a given azimuth
     * @param [in] az - azimuth index
     * @param [in,out] slice - Matrix2D object containing the slice extracted
     * @param [in] missing_value
     */
    void read_vertical_slice(unsigned az, Matrix2D<T>& slice, double missing_value) const
    {
        unsigned size_z = std::max(this->size(), (size_t)slice.rows());
        for (unsigned el = 0; el < size_z; ++el)
            this->scan(el).read_beam(az, slice.row_ptr(el), slice.cols(), missing_value);
    }

    /// Compute Volume statistics
    void compute_stats(VolumeStats& stats) const
    {
        stats.count_zeros.resize(this->size());
        stats.count_ones.resize(this->size());
        stats.count_others.resize(this->size());
        stats.sum_others.resize(this->size());

        for (unsigned iel = 0; iel < this->size(); ++iel)
        {
            stats.count_zeros[iel] = 0;
            stats.count_ones[iel] = 0;
            stats.count_others[iel] = 0;
            stats.sum_others[iel] = 0;

            for (unsigned iaz = 0; iaz < this->scan(iel).beam_count; ++iaz)
            {
                for (size_t i = 0; i < this->scan(iel).beam_size; ++i)
                {
                    int val = algo::DBZ::DBtoBYTE(this->scan(iel).get(iaz, i));
                    switch (val)
                    {
                        case 0: stats.count_zeros[iel]++; break;
                        case 1: stats.count_ones[iel]++; break;
                        default:
                                stats.count_others[iel]++;
                                stats.sum_others[iel] += val;
                                break;
                    }
                }
            }
        }
    }

    /**
     * Append a scan to this volume.
     *
     * It is required that scans are added in increasing elevation order,
     * because higher scan indices need to correspond to higher elevation
     * angles.
     * @param [in] beam_size - beam_size dimension
     * @param [in] elevation - PolarScan elevation [degrees]
     * @param [in] cell_size - cell_size value [m]
     * @return reference to the PolaScan added
     */
    PolarScan<T>& append_scan(unsigned beam_size, double elevation, double cell_size)
    {
        return volume::Scans<T>::append_scan(beam_count, beam_size, elevation, cell_size);
    }

    /**
     * Create or reuse a scan at position idx, with the given beam size
     * @param [in] idx - PolarScan Index 
     * @param [in] beam_size - beam_size dimension
     * @param [in] elevation - PolarScan elevation [degrees]
     * @param [in] cell_size - cell_size value [m]
     * @return reference to the PolaScan added
     */
    PolarScan<T>& make_scan(unsigned idx, unsigned beam_size, double elevation, double cell_size)
    {
        return volume::Scans<T>::make_scan(idx, beam_count, beam_size, elevation, cell_size);
    }

    /**
     * *= operator defined
     * @param [in] coefficient  - multiplicative constant
     */
    Volume& operator*=(const T coefficient)
    {
        for(unsigned el=0;el<this->size();el++)
        {
            this->scan(el)*=coefficient;
        }
        return *this;
    }

    /**
     * += operator defined
     * @param [in] addend - additive constant
     */
    Volume& operator+=(Volume& addend)
    {
        for(unsigned el=0;el<this->size();el++)
        {
            this->scan(el)+=addend[el];
        }
        return *this;
    }

protected:
    void resize_elev_fin();

private:
    Volume& operator=(Volume&);
    Volume(const Volume&);
};

/*
 * Tell the compiler that the implementation of these templates is already
 * explicitly instantiated in volume.cpp, and it should not spend time
 * reinstantiating it every time the template is used.
 */
extern template class PolarScan<double>;
extern template class Volume<double>;
namespace volume {
extern template class Scans<double>;
}

}

#endif
