#ifndef ARCHIVIATORE_VOLUME_CLASS_H
#define ARCHIVIATORE_VOLUME_CLASS_H

#include <logging.h>
#include <matrix.h>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

// TODO: prima o poi arriviamo a far senza di questi define
#define NUM_AZ_X_PPI 400

#define MISSING_DB (-20.)
#define MINVAL_DB (80./255.-20.)
#define MAXVAL_DB (60.)

namespace cumbac {

inline double BYTEtoDB(unsigned char z)
{
    return (z*80./255.-20.);
}

inline unsigned char DBtoBYTE(double dB)
{
    int byt = round((dB+20.)*255./80.);
    if (byt >= 0 && byt <= 255)
        return ((unsigned char)byt);
    else if (byt < 0)
        return 0;
    else
        return 255;
}


template<typename T>
class PolarScan : public Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
public:
    /// Count of beams in this scan
    const unsigned beam_count;
    /// Number of samples in each beam
    const unsigned beam_size;
    /**
     * Nominal elevation of this PolarScan, which may be different from the
     * effective elevation of each single beam
     */
    double elevation;

    PolarScan(unsigned beam_size, const T& default_value = BYTEtoDB(1))
        : Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(
                PolarScan::Constant(NUM_AZ_X_PPI, beam_size, default_value)), beam_count(NUM_AZ_X_PPI), beam_size(beam_size), elevation(0)
    {
    }

    ~PolarScan()
    {
    }

    /// Get a beam value
    T get(unsigned az, unsigned beam) const
    {
        return (*this)(az, beam);
    }

    /// Set a beam value
    void set(unsigned az, unsigned beam, T val)
    {
        (*this)(az, beam) = val;
    }

    /**
     * Riempie un array di float con i dati del raggio convertiti in DB
     * Se l'array è più lungo del raggio, setta gli elementi extra a missing.
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
};

struct VolumeStats
{
    std::vector<unsigned> count_zeros;
    std::vector<unsigned> count_ones;
    std::vector<unsigned> count_others;
    std::vector<unsigned> sum_others;

    void print(FILE* out);
};

template<typename T>
class Volume : protected std::vector<PolarScan<T>*>
{
public:
    using std::vector<PolarScan<T>*>::size;
    typedef typename std::vector<PolarScan<T>*>::iterator iterator;
    typedef typename std::vector<PolarScan<T>*>::const_iterator const_iterator;

    // Access a polar scan
    PolarScan<T>& scan(unsigned idx) { return *(*this)[idx]; }
    const PolarScan<T>& scan(unsigned idx) const { return *(*this)[idx]; }

    Volume()
    {
    }

    template<typename OT>
    Volume(const Volume<OT>& v, const T& default_value)
    {
        this->resize(v.size(), 0);
        for (unsigned i = 0; i < v.size(); ++i)
        {
            (*this)[i] = new PolarScan<T>(v.scan(i).beam_size, default_value);
            (*this)[i]->elevation = v.scan(i).elevation;
        }
    }

    ~Volume()
    {
        for (auto i: *this)
            if (i) delete i;
    }

    /// Return the maximum beam count in all PolarScans
    const unsigned max_beam_count() const
    {
        unsigned res = 0;
        for (size_t i = 0; i < size(); ++i)
            res = std::max(res, (*this)[i]->beam_count);
        return res;
    }

    /// Return the maximum beam size in all PolarScans
    const unsigned max_beam_size() const
    {
        unsigned res = 0;
        for (size_t i = 0; i < this->size(); ++i)
            res = std::max(res, (*this)[i]->beam_size);
        return res;
    }


    double elevation_min() const
    {
        return this->front()->elevation;
    }

    double elevation_max() const
    {
        return this->back()->elevation;
    }

    /**
     * Fill a matrix elevations x beam_size with the vertical slice at a given
     * azimuth
     */
    void read_vertical_slice(unsigned az, Matrix2D<T>& slice, double missing_value) const
    {
        unsigned size_z = std::max(size(), (size_t)slice.rows());
        for (unsigned el = 0; el < size_z; ++el)
            scan(el).read_beam(az, slice.row(el), slice.cols(), missing_value);
    }

    void compute_stats(VolumeStats& stats) const
    {
        stats.count_zeros.resize(this->size());
        stats.count_ones.resize(this->size());
        stats.count_others.resize(this->size());
        stats.sum_others.resize(this->size());

        for (int iel = 0; iel < this->size(); ++iel)
        {
            stats.count_zeros[iel] = 0;
            stats.count_ones[iel] = 0;
            stats.count_others[iel] = 0;
            stats.sum_others[iel] = 0;

            for (unsigned iaz = 0; iaz < scan(iel).beam_count; ++iaz)
            {
                for (size_t i = 0; i < scan(iel).beam_size; ++i)
                {
                    int val = DBtoBYTE(scan(iel).get(iaz, i));
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

    // Create or reuse a scan at position idx, with the given beam size
    PolarScan<T>& make_scan(unsigned idx, unsigned beam_size, double elevation)
    {
        if (idx < this->size())
        {
            if (beam_size != (*this)[idx]->beam_size)
            {
                LOG_CATEGORY("radar.io");
                LOG_ERROR("make_scan(idx=%u, beam_size=%u) called, but the scan already existed with beam_size=%u", idx, beam_size, (*this)[idx]->beam_size);
                throw std::runtime_error("beam_size mismatch");
            }
        } else {
            // If some elevation has been skipped, fill in the gap
            if (idx > this->size())
            {
                if (this->empty())
                    this->push_back(new PolarScan<T>(beam_size));
                while (this->size() < idx)
                    this->push_back(new PolarScan<T>(this->back()->beam_size));
            }

            // Add the new polar scan
            this->push_back(new PolarScan<T>(beam_size));
            this->back()->elevation = elevation;
        }

        // Return it
        return *(*this)[idx];
    }

protected:
    void resize_elev_fin();

private:
    Volume(const Volume&);
    Volume& operator=(const Volume&);
};


}

#endif
