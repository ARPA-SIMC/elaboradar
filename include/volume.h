#ifndef ARCHIVIATORE_VOLUME_CLASS_H
#define ARCHIVIATORE_VOLUME_CLASS_H

#include <string>
#include <vector>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <matrix.h>

// TODO: prima o poi arriviamo a far senza di questi define
#define NUM_AZ_X_PPI 400

namespace H5 {
struct H5File;
}

namespace cumbac {
struct Site;

namespace volume {

template<typename T> struct WithLoadInfo;

class LoadInfo
{
private:
    /**
     * Reference counter, to share instances across many volumes
     */
    unsigned _refcount;

protected:
    LoadInfo() : _refcount(0) {}
    virtual ~LoadInfo() {}

    // Increment the reference count
    void _ref() { ++_refcount; }

    /**
     * Decrement the reference count, and return true if it became zero and
     * this structure should be deallocated
     */
    bool _unref() { return --_refcount == 0; }

    template <typename T> friend struct WithLoadInfo;
};

template<class LI>
struct WithLoadInfo
{
    LI* _load_info;

    WithLoadInfo() : _load_info(0) {}
    WithLoadInfo(const WithLoadInfo<LI>& o)
        : _load_info(o._load_info ? o._load_info->_ref() : 0)
    {
    }
    ~WithLoadInfo()
    {
        if (_load_info && _load_info->_unref())
            delete _load_info;
    }
    WithLoadInfo<LI> operator=(const WithLoadInfo& o)
    {
        if (_load_info && _load_info != o._load_info && _load_info->_unref())
            delete _load_info;
        _load_info = o._load_info;
        _load_info->_ref();
        return *this;
    }

    LI* get_load_info() const { return _load_info; }

    LI& obtain_load_info()
    {
        if (!_load_info)
        {
            _load_info = new LI();
            _load_info->_ref();
        }
        return *_load_info;
    }

    const LI& load_info() const
    {
        if (!_load_info) throw std::runtime_error("this object has no loading information");
        return *_load_info;
    }

    LI* add_load_info(LI* info)
    {
        if (_load_info && _load_info != info && _load_info->_unref())
            delete _load_info;
        _load_info = info;
        _load_info->_ref();
        return _load_info;
    }
};

struct LoadLogEntry
{
    double theta;
    double alpha;

    LoadLogEntry(double theta, double alpha)
        : theta(theta), alpha(alpha)
    {
    }

    bool operator==(const LoadLogEntry& e) const
    {
        return theta == e.theta && alpha == e.alpha;
    }
};

struct LoadLog : public std::vector<LoadLogEntry>
{
    void log(double theta, double alpha)
    {
        push_back(LoadLogEntry(theta, alpha));
    }
    void print(FILE* out) const;
};

struct BeamInfo
{
    /// Real beam elevation in degrees
    LoadLog load_log;
    double elevation;
};

}

struct PolarScanLoadInfo : public volume::LoadInfo
{
    std::vector<volume::BeamInfo> beam_info;

    PolarScanLoadInfo(unsigned beam_count)
    {
        beam_info.resize(beam_count);
    }

    /// Return the number of beams that have been filled with data while loading
    unsigned count_rays_filled() const
    {
        unsigned count = 0;
        for (std::vector<volume::BeamInfo>::const_iterator i = beam_info.begin(); i != beam_info.end(); ++i)
            if (!i->load_log.empty())
                ++count;
        return count;
    }

    /// Return the load log for the given beam
    const volume::LoadLog& get_beam_load_log(unsigned az) const
    {
        return beam_info[az].load_log;
    }

    inline double get_elevation(unsigned az) const
    {
        return beam_info[az].elevation;
    }

    inline double get_elevation_rad(unsigned az) const
    {
        return beam_info[az].elevation * M_PI / 180.;
    }
};

struct VolumeLoadInfo : public volume::LoadInfo
{
    std::string filename;
    bool declutter_rsp; // ?

    VolumeLoadInfo()
        : declutter_rsp(false)
    {
    }
};

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
class PolarScan : public Matrix2D<T>, public volume::WithLoadInfo<PolarScanLoadInfo>
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

    PolarScan(unsigned beam_size)
        : Matrix2D<T>(beam_size, NUM_AZ_X_PPI, BYTEtoDB(1)), beam_count(NUM_AZ_X_PPI), beam_size(beam_size), elevation(0)
    {
    }

    ~PolarScan()
    {
    }

    /// Get a raw value in a beam
    unsigned char get_raw(unsigned az, unsigned beam) const
    {
        return DBtoBYTE(get_db(az, beam));
    }

    /// Get a beam value in DB
    double get_db(unsigned az, unsigned beam) const
    {
        return (*this)[az][beam];
    }

    /// Set a raw value in a beam
    void set_raw(unsigned az, unsigned beam, unsigned char val)
    {
        (*this)[az][beam] = BYTEtoDB(val);
    }

    /// Set a beam value in DB
    void set_db(unsigned az, unsigned beam, double val)
    {
        (*this)[az][beam] = val;
    }

    /**
     * Riempie un array di float con i dati del raggio convertiti in DB
     * Se l'array è più lungo del raggio, setta gli elementi extra a missing.
     */
    void read_beam_db(unsigned az, float* out, unsigned out_size, float missing=0) const
    {
        using namespace std;

        // Prima riempio il minimo tra ray.size() e out_size
        size_t set_count = min(beam_size, out_size);

        for (unsigned i = 0; i < set_count; ++i)
            out[i] = get_db(az, i);

        for (unsigned i = set_count; i < out_size; ++i)
            out[i] = missing;
    }


    void fill_beam(int el_num, double theta, double alpha, unsigned size, const double* data);


protected:
    void merge_beam(int el_num, int az_num, double theta, double alpha, unsigned size, const double* dati);
};

struct VolumeStats
{
    std::vector<unsigned> count_zeros;
    std::vector<unsigned> count_ones;
    std::vector<unsigned> count_others;
    std::vector<unsigned> sum_others;

    void print(FILE* out);
};

struct VolumeLoadOptions
{
    const Site& site;
    bool medium;
    bool clean;
    std::vector<double> elev_array;
    /**
     * If this is greather than zero, truncate each beam to this number of
     * samples
     */
    unsigned max_bin;

    VolumeLoadOptions(const Site& site, bool medium=false, bool clean=false, unsigned max_bin=0);

    /**
     * Compute the vol_pol index of an elevation angle
     * @returns -1 if no suitable index was found, else the index
     */
    int elevation_index(double elevation) const;
};

template<typename T>
class Volume : public volume::WithLoadInfo<VolumeLoadInfo>
{
protected:
    // Dato di base volume polare
    std::vector<PolarScan<T>*> scans;

    // Create or reuse a scan at position idx, with the given beam size
    PolarScan<T>& make_scan(const VolumeLoadOptions& opts, unsigned idx, unsigned beam_size);

public:
    // Acquisition date
    time_t acq_date;
    // Length of a beam cell in meters
    double size_cell;
    // How many elevations are in this volume
    unsigned int NEL;

    // Access a polar scan
    PolarScan<T>& scan(unsigned idx) { return *scans[idx]; }
    const PolarScan<T>& scan(unsigned idx) const { return *scans[idx]; }

    // elevazione finale in coordinate azimut range
    std::vector<unsigned char> elev_fin[NUM_AZ_X_PPI];

    Volume()
        : acq_date(0), size_cell(0), NEL(0)
    {
    }

    ~Volume()
    {
        for (typename std::vector<PolarScan<T>*>::iterator i = scans.begin(); i != scans.end(); ++i)
            if (*i) delete *i;
    }

    /// Return the maximum beam count in all PolarScans
    const unsigned max_beam_count() const
    {
        unsigned res = 0;
        for (size_t i = 0; i < scans.size(); ++i)
            res = std::max(res, scans[i]->beam_count);
        return res;
    }

    /// Return the maximum beam size in all PolarScans
    const unsigned max_beam_size() const
    {
        unsigned res = 0;
        for (size_t i = 0; i < scans.size(); ++i)
            res = std::max(res, scans[i]->beam_size);
        return res;
    }


    double elevation_min() const
    {
        return scans.front()->elevation;
    }

    double elevation_max() const
    {
        return scans.back()->elevation;
    }

    inline double elevation_rad_at_elev_preci(unsigned az_idx, unsigned ray_idx) const
    {
        return scan(elev_fin[az_idx][ray_idx]).load_info().get_elevation_rad(az_idx);
    }

    inline double elevation_at_elev_preci(unsigned az_idx, unsigned ray_idx) const
    {
        return scan(elev_fin[az_idx][ray_idx]).load_info().get_elevation(az_idx);
    }

    inline unsigned char sample_at_elev_preci(unsigned az_idx, unsigned ray_idx) const
    {
        const PolarScan<T>& s = scan(elev_fin[az_idx][ray_idx]);
        if (ray_idx < s.beam_size)
            return s.get_raw(az_idx, ray_idx);
        else
            // If we are reading out of bounds, return 1 (the missing value)
            return 1;
    }

    void read_sp20(const char* nome_file, const VolumeLoadOptions& options);
    void read_odim(const char* nome_file, const VolumeLoadOptions& options);

    void compute_stats(VolumeStats& stats) const;

    void write_info_to_debug_file(H5::H5File out);

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
