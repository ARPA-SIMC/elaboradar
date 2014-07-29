#ifndef ARCHIVIATORE_VOLUME_CLASS_H
#define ARCHIVIATORE_VOLUME_CLASS_H

#include <logging.h>
#include <matrix.h>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <memory>
#include <Eigen/Core>

// TODO: prima o poi arriviamo a far senza di questi define
#define NUM_AZ_X_PPI 400

#define MISSING_DB (-20.)
#define MINVAL_DB (80./255.-20.)
#define MAXVAL_DB (60.)

namespace elaboradar {

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
    unsigned beam_count;
    /// Number of samples in each beam
    unsigned beam_size;
    /**
     * Nominal elevation of this PolarScan, which may be different from the
     * effective elevation of each single beam
     */
    double elevation;
    /// Vector of actual elevations for each beam
    Eigen::VectorXd elevations_real;
    /// Size of a beam cell in meters
    double cell_size;
    T nodata = 0;    // Value used as 'no data' value
    T undetect = 0;  // Minimum amount that can be measured
    T gain = 1;
    T offset = 0;

    PolarScan(unsigned beam_count, unsigned beam_size, const T& default_value = BYTEtoDB(1))
        : Matrix2D<T>(PolarScan::Constant(beam_count, beam_size, default_value)),
          beam_count(beam_count), beam_size(beam_size), elevation(0), elevations_real(beam_count)
    {
    }

    PolarScan(const PolarScan& s) = default;
#if 0
    // PolarScan(PolarScan&& s) = default;
    PolarScan(const PolarScan& s)
        : Matrix2D<T>(PolarScan::Constant(s.beam_count, s.beam_size, s.nodata)),
          beam_count(s.beam_count), beam_size(s.beam_size),
          elevation(s.elevation), elevations_real(s.elevations_real),
          nodata(s.nodata), undetect(s.undetect), gain(s.gain), offset(s.offset)
    {
    }
#endif

    template<typename OT>
    PolarScan(const PolarScan<OT>& s, const T& default_value)
        : Matrix2D<T>(PolarScan::Constant(s.beam_count, s.beam_size, default_value)),
          beam_count(s.beam_count), beam_size(s.beam_size),
          elevation(s.elevation), elevations_real(s.elevations_real),
          nodata(default_value)
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

namespace volume {

struct LoadInfo
{
    std::string filename;
    // Acquisition date
    time_t acq_date;
    bool declutter_rsp; // ?

    LoadInfo()
        : declutter_rsp(false)
    {
    }
};

/**
 * Volume of polarscans which can have a different beam count for each elevation
 */
template<typename T>
class Scans : public std::vector<PolarScan<T>>
{
public:
    typedef T Scalar;
    std::string quantity;
    std::string units;
    std::shared_ptr<LoadInfo> load_info;

    Scans() = default;

    template<typename OT>
    Scans(const Scans<OT>& v, const T& default_value)
    {
        this->quantity = v.quantity;
        this->units = v.units;
        this->load_info = v.load_info;
        this->reserve(v.size());

        for (const auto& src_scan : v)
            this->push_back(PolarScan<T>(src_scan, default_value));
    }

    // Create or reuse a scan at position idx, with the given beam size
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
};

}

/**
 * Homogeneous volume with a common beam count for all PolarScans
 */
template<typename T>
class Volume : public volume::Scans<T>
{
public:
    const unsigned beam_count;

    // Access a polar scan
    PolarScan<T>& scan(unsigned idx) { return (*this)[idx]; }
    const PolarScan<T>& scan(unsigned idx) const { return (*this)[idx]; }

    Volume(unsigned beam_count=NUM_AZ_X_PPI)
        : beam_count(beam_count)
    {
    }

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


    double elevation_min() const
    {
        return this->front().elevation;
    }

    double elevation_max() const
    {
        return this->back().elevation;
    }

    /**
     * Fill a matrix elevations x beam_size with the vertical slice at a given
     * azimuth
     */
    void read_vertical_slice(unsigned az, Matrix2D<T>& slice, double missing_value) const
    {
        unsigned size_z = std::max(this->size(), (size_t)slice.rows());
        for (unsigned el = 0; el < size_z; ++el)
            scan(el).read_beam(az, slice.row_ptr(el), slice.cols(), missing_value);
    }

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
    PolarScan<T>& make_scan(unsigned idx, unsigned beam_size, double elevation, double cell_size)
    {
        return volume::Scans<T>::make_scan(idx, beam_count, beam_size, elevation, cell_size);
    }

    void filter(Volume<T>& raw, double filter_range)
    {
	unsigned window_size;
	this->quantity=raw.quantity;
	this->clear();
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*filter_range/raw.scan(i).cell_size);
		this->make_filter_scan_range(raw.scan(i),window_size);
	}
    }

    void textureSD(Volume<T>& raw, double filter_range)
    {
	Volume<T> filtered;
	filtered.filter(raw,filter_range);
	unsigned window_size;
	this->quantity=raw.quantity;
	this->clear();
	for(unsigned i=0;i<raw.size();i++)
	{
		window_size=1+2*std::floor(0.5*filter_range/raw.scan(i).cell_size);
		this->make_rms_scan_range(raw.scan(i),filtered.scan(i),window_size);
	}
    }

    void lin2dB(Volume<T>& lin)
    {
	//this->quantity=lin.quantity.lin2dB(); // TODO: not yet implemented
	this->clear();
	for(unsigned i=0;i<lin.size();i++)
	{
		this->push_back(PolarScan<T>(lin.scan(i).beam_count,lin.scan(i).beam_size,0.));
		this->back().elevation = lin.scan(i).elevation;
		this->back().cell_size = lin.scan(i).cell_size;
		this->back().block(0,0,lin.scan(i).beam_count,lin.scan(i).beam_size)=lin.scan(i).log10();		
		this->back().array()*=10.;
	}
    }

    void dB2lin(Volume<T>& DB)
    {
	//this->quantity=DB.quantity.dB2lin(); // TODO: not yet implemented
	this->clear();
	for(unsigned i=0;i<DB.size();i++)
	{
		this->push_back(PolarScan<T>(DB.scan(i).beam_count,DB.scan(i).beam_size,0.));
		this->back().elevation = DB.scan(i).elevation;
		this->back().cell_size = DB.scan(i).cell_size;
		this->back().block(0,0,DB.scan(i).beam_count,DB.scan(i).beam_size)=DB.scan(i).array()*0.1;
		this->back().block(0,0,DB.scan(i).beam_count,DB.scan(i).beam_size)=this->back().exp10();
	}
    }

    Volume& operator*=(const T coefficient)
    {
	for(unsigned el=0;el<this->size();el++)
	{
		this->scan(el)*=coefficient;
	}
	return *this;
    }

    Volume& operator=(const Volume&){return *this;};
    Volume(const Volume&);

protected:
    void resize_elev_fin();

private:
   
    

    /**
     * Fill this volume by filtering through a mobile window the data from an
     * existing volume
     */
    void make_filter_scan_range(PolarScan<T>& raw, unsigned win)
    {
	unsigned half_win=0.5*(win-1);
	this->push_back(PolarScan<T>(raw.beam_count,raw.beam_size,0.));
	this->back().elevation = raw.elevation;
	this->back().cell_size = raw.cell_size;
	Matrix2D<unsigned> counter(Matrix2D<unsigned>::Constant(raw.beam_count,raw.beam_size,0));
	T value;
	for(unsigned i=0;i<raw.rows();i++)
	{
		for(unsigned j=0;j<half_win;j++)
		{
			value=raw(i,j);
			if((value != raw.undetect) && (value != raw.nodata))
			{
				this->back().block(i,0,1,half_win+1+j).array()+=value;
				counter.block(i,0,1,half_win+1+j).array()+=1;
			}
			value=raw(i,raw.beam_size-half_win+j);
			if((value != raw.undetect)&&(value != raw.nodata))
			{
				this->back().block(i,raw.beam_size-win+1+j,1,win-1-j).array()+=value;
				counter.block(i,raw.beam_size-win+1+j,1,win-1-j).array()+=1;
			}
		}
		for(unsigned j=half_win;j<(raw.beam_size-half_win);j++)
		{
			value=raw(i,j);
			if((value!= raw.undetect) && (value!=raw.nodata))
			{
				this->back().block(i,j-half_win,1,win).array()+=value;
				counter.block(i,j-half_win,1,win).array()+=1;
			}
		}
		for(unsigned j=0;j<raw.beam_size;j++)
		{
			if(counter(i,j)) this->back().set(i,j,(this->back().get(i,j))/counter(i,j));
			else this->back().set(i,j, raw.undetect);
		}
	}
    }

    void make_rms_scan_range(PolarScan<T>& raw, PolarScan<T>& filtered, unsigned win)
    {
	unsigned half_win=0.5*(win-1);
	this->push_back(PolarScan<T>(raw.beam_count,raw.beam_size,0.));
	this->back().elevation = raw.elevation;
	this->back().cell_size = raw.cell_size;
	Matrix2D<unsigned> counter(Matrix2D<unsigned>::Constant(beam_count,raw.beam_size,0));
	T value;
	for(unsigned i=0;i<raw.rows();i++)
	{
		for(unsigned j=0;j<half_win;j++)
		{
			if((raw(i,j)!=raw.undetect)&&(raw(i,j)!=raw.nodata))
			{
				value=(raw(i,j)-filtered(i,j))*(raw(i,j)-filtered(i,j));
				this->back().block(i,0,1,half_win+1+j).array()+=value;
				counter.block(i,0,1,half_win+1+j).array()+=1;
			}
			if((raw(i,raw.beam_size-half_win+j)!=raw.undetect)&&(raw(i,raw.beam_size-half_win+j)!=raw.nodata))
			{
				value=(raw(i,raw.beam_size-half_win+j)-filtered(i,raw.beam_size-half_win+j));
				this->back().block(i,raw.beam_size-win+1+j,1,win-1-j).array()+=value;
				counter.block(i,raw.beam_size-win+1+j,1,win-1-j).array()+=1;
			}
		}
		for(unsigned j=half_win;j<(raw.beam_size-half_win);j++)
		{
			if((raw(i,j)!=raw.undetect)&&(raw(i,j)!=raw.nodata))
			{
				value=(raw(i,j)-filtered(i,j))*(raw(i,j)-filtered(i,j));
				this->back().block(i,j-half_win,1,win).array()+=value;
				counter.block(i,j-half_win,1,win).array()+=1;
			}
		}
		for(unsigned j=0;j<raw.beam_size;j++)
		{
			if(counter(i,j)) this->back().set(i,j,std::sqrt(this->back().get(i,j)/counter(i,j)));
			else this->back().set(i,j,raw.undetect);
		}
	}
    }

    
};


}

#endif
