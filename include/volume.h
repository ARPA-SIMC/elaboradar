#ifndef ARCHIVIATORE_VOLUME_CLASS_H
#define ARCHIVIATORE_VOLUME_CLASS_H

#include <logging.h>
#include <matrix.h>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>

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
class PolarScan : public Matrix2D<T>
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
        : Matrix2D<T>(beam_size, NUM_AZ_X_PPI, default_value), beam_count(NUM_AZ_X_PPI), beam_size(beam_size), elevation(0)
    {
    }

    ~PolarScan()
    {
    }

    /// Get a beam value
    T get(unsigned az, unsigned beam) const
    {
        return (*this)[az][beam];
    }

    /// Set a beam value
    void set(unsigned az, unsigned beam, T val)
    {
        (*this)[az][beam] = val;
    }

    /**
     * Riempie un array di float con i dati del raggio convertiti in DB
     * Se l'array è più lungo del raggio, setta gli elementi extra a missing.
     */
    void read_beam_db(unsigned az, double* out, unsigned out_size, float missing=0) const
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
    // How many elevations are in this volume
    unsigned int NEL;

    using std::vector<PolarScan<T>*>::size;
    typedef typename std::vector<PolarScan<T>*>::iterator iterator;
    typedef typename std::vector<PolarScan<T>*>::const_iterator const_iterator;

    // Access a polar scan
    PolarScan<T>& scan(unsigned idx) { return *(*this)[idx]; }
    const PolarScan<T>& scan(unsigned idx) const { return *(*this)[idx]; }

    Volume()
        : NEL(0)
    {
    }

    template<typename OT>
    Volume(const Volume<OT>& v, const T& default_value)
        : NEL(v.NEL)
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
        for (iterator i = this->begin(); i != this->end(); ++i)
            if (*i) delete *i;
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
        // Enlarge the scans vector if needed
        if (idx >= this->size())
        {
            this->resize(idx + 1, 0);
            this->NEL = idx + 1;
        }

        // Create the PolarScan if needed
        if (!(*this)[idx])
        {
            (*this)[idx] = new PolarScan<T>(beam_size);
            (*this)[idx]->elevation = elevation;
        }
        else if (beam_size != (*this)[idx]->beam_size)
        {
            LOG_CATEGORY("radar.io");
            LOG_ERROR("make_scan(idx=%u, beam_size=%u) called, but the scan already existed with beam_size=%u", idx, beam_size, (*this)[idx]->beam_size);
            throw std::runtime_error("beam_size mismatch");
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


template<typename T>
struct ArrayStats
{
    bool first;
    T min;
    T max;
    double avg;
    unsigned count_zeros;
    unsigned count_ones;

    ArrayStats() : first(true), count_zeros(0), count_ones(0) {}

    void count_sample(const T& sample, unsigned item_count)
    {
        if (sample == 0)
            ++count_zeros;
        else if (sample == 1)
            ++count_ones;

        if (first)
        {
            min = sample;
            max = sample;
            avg = (double)min / item_count;
            first = false;
        }
        else
        {
            if (sample < min)
                min = sample;
            if (sample > max)
                max = sample;
            avg += (double)sample / item_count;
        }
    }
};

/**
 * Store one element for each sample in a volume
 */
template<typename T>
class VolumeInfo
{
protected:
    const unsigned sz_el;
    const unsigned sz_az;
    const unsigned sz_beam;
    T* data;

public:
    template<typename BIN>
    VolumeInfo(const Volume<BIN>& vol)
        : sz_el(vol.NEL), sz_az(vol.max_beam_count()), sz_beam(vol.max_beam_size()),
          data(new T[sz_el * sz_az * sz_beam])
    {
    }

    ~VolumeInfo()
    {
        delete[] data;
    }

    // Fill all the volume info with the given value
    void init(const T& val)
    {
        for (unsigned i = 0; i < sz_el * sz_az * sz_beam; ++i)
            data[i] = val;
    }

    const T& get(unsigned el, unsigned az, unsigned beam) const
    {
        return data[el * (sz_az * sz_beam) + az * sz_beam + beam];
    }

    void set(unsigned el, unsigned az, unsigned beam, const T& val)
    {
        data[el * (sz_az * sz_beam) + az * sz_beam + beam] = val;
    }

    void fill_array_stats(ArrayStats<T>& stats) const
    {
        for (unsigned i = 0; i < sz_el * sz_az * sz_beam; ++i)
            stats.count_sample(data[i], sz_el * sz_az * sz_beam);
    }
};

}

#endif
