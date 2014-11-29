#ifndef ARCHIVIATORE_TEST_UTILS_CLASS_H
#define ARCHIVIATORE_TEST_UTILS_CLASS_H

#include <wibble/tests.h>
#include <elaboradar/volume.h>

namespace elaboradar {
struct Cart;
struct CartLowris;
struct CUM_BAC;
}

namespace testradar {

template<typename T>
struct ArrayStats
{
    bool all_missing = true;
    T min = 0;
    T max = 0;
    double sum = 0;
    unsigned count_missing = 0;
    unsigned count_values = 0;

    ArrayStats() {}

    double avg() const { return count_values ? sum / count_values : min; }

    void count_sample(const T& sample)
    {
        if (all_missing)
        {
            min = sample;
            max = sample;
            all_missing = false;
        }
        else
        {
            if (sample < min)
                min = sample;
            if (sample > max)
                max = sample;
        }
        sum += (double)sample;
        ++count_values;
    }

    void count_sample(const T& missing, const T& sample)
    {
        if (sample == missing)
            ++count_missing;
        else
            count_sample(sample);
    }


    void fill(const T* arr, unsigned size)
    {
        for (unsigned i = 0; i < size; ++i)
            this->count_sample(arr[i]);
    }

    void fill(const elaboradar::Matrix2D<T>& arr)
    {
        this->fill(arr.data(), arr.size());
    }

    void fill(const elaboradar::Volume<T>& vol)
    {
        for (unsigned i = 0; i < vol.size(); ++i)
            this->fill(vol.scan(i));
    }


    void fill(const T& missing, const T* arr, unsigned size)
    {
        for (unsigned i = 0; i < size; ++i)
            this->count_sample(missing, arr[i]);
    }

    void fill(const T& missing, const elaboradar::Matrix2D<T>& arr)
    {
        this->fill(missing, arr.data(), arr.size());
    }

    void fill(const T& missing, const elaboradar::Volume<T>& vol)
    {
        for (unsigned i = 0; i < vol.size(); ++i)
            this->fill(missing, vol.scan(i));
    }


    void print()
    {
//        fprintf(stderr, "min %f max %f avg  %f, zeros: %u, ones: %u\n",
//                (double)this->min, (double)this->max, this->avg(),
//                this->count_zeros, this->count_ones);
        fprintf(stderr, "min %f max %f \n",(double)this->min, (double)this->max);
    }
};

template<typename T> inline T to_num(const T& val) { return val; }
inline double to_num(const double& val) { return round(val * 100.0) / 100.0; }
inline float to_num(const float& val) { return round(val * 100.0) / 100.0; }
inline unsigned to_num(const unsigned char& val) { return val; }
inline int to_num(const char& val) { return val; }

template<typename T> inline bool approx_equals(const T& v1, const T& v2) { return v1 == v2; }
inline bool approx_equals(const double& v1, const double& v2) { return round(v1 * 100) == round(v2 * 100); }
inline bool approx_equals(const float& v1, const float& v2) { return roundf(v1 * 100) == roundf(v2 * 100); }

template<typename DATA>
struct TestStatsEqual
{
    typedef typename DATA::Scalar Scalar;
    const DATA& matrix;
    bool has_missing = false;
    Scalar missing;
    unsigned count_missing = 0;
    Scalar min;
    double avg;
    Scalar max;

    TestStatsEqual(const DATA& actual, Scalar min, double avg, Scalar max)
        : matrix(actual), min(min), avg(avg), max(max)
    {
    }

    TestStatsEqual(const DATA& actual, Scalar missing, unsigned count_missing, Scalar min, double avg, Scalar max)
        : matrix(actual), has_missing(true), missing(missing), count_missing(count_missing), min(min), avg(avg), max(max)
    {
    }

    void check(WIBBLE_TEST_LOCPRM) const
    {
        using namespace wibble::tests;
        using namespace std;

        ArrayStats<Scalar> stats;
        bool failed = false;
        if (has_missing)
        {
            stats.fill(missing, matrix);
            if (stats.count_missing != count_missing)
                failed = true;
        } else
            stats.fill(matrix);
        if (!approx_equals(stats.min, min)) failed = true;
        if (!approx_equals(stats.max, max)) failed = true;
        if (!approx_equals(stats.avg(), avg)) failed = true;

        if (failed)
        {
            std::stringstream ss;
            ss << "stats (";
            if (has_missing)
                ss << "missing: " << stats.count_missing << " ";
            ss << "min: " << to_num(stats.min)
               << " avg: " << to_num(stats.avg())
               << " max: " << to_num(stats.max)
               << ") differ from expected (";
            if (has_missing)
                ss << "missing: " << count_missing << " ";
            ss << "min: " << to_num(min)
               << " avg: " << to_num(avg)
               << " max: " << to_num(max)
               << ")";
            wibble_test_location.fail_test(ss.str());
        }
    }
};

template<typename T>
struct ActualMatrix2D : public wibble::tests::Actual<const elaboradar::Matrix2D<T>&>
{
    using wibble::tests::Actual<const elaboradar::Matrix2D<T>&>::Actual;

    template<typename... args>
    TestStatsEqual<elaboradar::Matrix2D<T>> statsEqual(args&&... params) const
    {
        return TestStatsEqual<elaboradar::Matrix2D<T>>(this->actual, params...);
    }
};

template<typename T>
struct ActualVolume : public wibble::tests::Actual<const elaboradar::Volume<T>&>
{
    using wibble::tests::Actual<const elaboradar::Volume<T>&>::Actual;

    template<typename... args>
    TestStatsEqual<elaboradar::Volume<T>> statsEqual(args&&... params) const
    {
        return TestStatsEqual<elaboradar::Volume<T>>(this->actual, params...);
    }
};

template<typename T>
inline ActualMatrix2D<T> actual(const elaboradar::Matrix2D<T>& actual) { return ActualMatrix2D<T>(actual); }

template<typename T>
inline ActualMatrix2D<T> actual(const elaboradar::PolarScan<T>& actual) { return ActualMatrix2D<T>(actual); }

template<typename T>
inline ActualMatrix2D<T> actual(const elaboradar::Image<T>& actual) { return ActualMatrix2D<T>(actual); }

template<typename T>
inline ActualVolume<T> actual(const elaboradar::Volume<T>& actual) { return ActualVolume<T>(actual); }

template<typename DATA>
void print_stats(const std::string& name, const DATA& data, const typename DATA::Scalar& missing, std::ostream& out)
{
    using namespace std;
    ArrayStats<typename DATA::Scalar> stats;
    stats.fill(missing, data);
    out << "wassert(actual(" << name << ").statsEqual"
        << "(" << to_num(missing)
        << ", " << stats.count_missing
        << ", " << to_num(stats.min)
        << ", " << to_num(stats.avg())
        << ", " << to_num(stats.max)
        << "));" << endl;
}

template<typename DATA>
void print_stats(const std::string& name, const DATA& data, std::ostream& out)
{
    using namespace std;
    ArrayStats<typename DATA::Scalar> stats;
    stats.fill(data);
    out << "wassert(actual(" << name << ").statsEqual"
        << "(" << to_num(stats.min)
        << ", " << to_num(stats.avg())
        << ", " << to_num(stats.max)
        << "));" << endl;
}

void print_stats(const std::string& name, const elaboradar::CUM_BAC& cb, std::ostream& out);
void print_stats(const std::string& name, const elaboradar::Cart& cart, std::ostream& out);
void print_stats(const std::string& name, const elaboradar::CartLowris& cart, std::ostream& out);

}

#endif
