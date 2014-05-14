#ifndef ARCHIVIATORE_TEST_UTILS_CLASS_H
#define ARCHIVIATORE_TEST_UTILS_CLASS_H

#include <wibble/tests.h>
#include <volume.h>
#include <iomanip>

namespace testradar {

template<typename T>
struct ArrayStats
{
    bool all_missing = false;
    T min = 0;
    T max = 0;
    double avg = 0;
    unsigned count_missing = 0;

    ArrayStats() {}

    void count_sample(const T& sample, unsigned item_count)
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
        avg += (double)sample / item_count;
    }

    void count_sample(const T& missing, const T& sample, unsigned item_count)
    {
        if (sample == missing)
            ++count_missing;
        else
            count_sample(sample, item_count);
    }


    void fill(const T* arr, unsigned size)
    {
        for (unsigned i = 0; i < size; ++i)
            this->count_sample(arr[i], size);
    }

    void fill(const cumbac::Matrix2D<T>& arr)
    {
        this->fill(arr.data(), arr.size());
    }

    void fill(const cumbac::Volume<T>& vol)
    {
        unsigned nsamples = 0;
        for (unsigned i = 0; i < vol.size(); ++i)
            nsamples += vol.scan(i).size();

        for (unsigned i = 0; i < vol.size(); ++i)
            for (size_t j = 0; j < vol.scan(i).size(); ++j)
                this->count_sample(vol.scan(i).data()[j], nsamples);
    }


    void fill(const T& missing, const T* arr, unsigned size)
    {
        for (unsigned i = 0; i < size; ++i)
            this->count_sample(missing, arr[i], size);
    }

    void fill(const T& missing, const cumbac::Matrix2D<T>& arr)
    {
        this->fill(missing, arr.data(), arr.size());
    }

    void fill(const T& missing, const cumbac::Volume<T>& vol)
    {
        unsigned nsamples = 0;
        for (unsigned i = 0; i < vol.size(); ++i)
            nsamples += vol.scan(i).size();

        for (unsigned i = 0; i < vol.size(); ++i)
            for (size_t j = 0; j < vol.scan(i).size(); ++j)
                this->count_sample(missing, vol.scan(i).data()[j], nsamples);
    }


    void print()
    {
        fprintf(stderr, "min %f max %f avg %f, zeros: %u, ones: %u\n",
                (double)this->min, (double)this->max, this->avg,
                this->count_zeros, this->count_ones);
    }
};

template<typename T>
int avg(const cumbac::Matrix2D<T>& m)
{
    const unsigned size = m.rows() * m.cols();
    double mean = 0;
    for (unsigned i = 0; i < size; ++i)
        mean += (double)m.data()[i] / size;
    return round(mean);
}

template<typename DATA, typename T>
struct TestStatsEqual
{
    const DATA& matrix;
    bool has_missing = false;
    unsigned count_missing = 0;
    T missing;
    T min;
    double avg;
    T max;

    TestStatsEqual(const DATA& actual, T min, double avg, T max)
        : matrix(actual), min(min), avg(avg), max(max)
    {
    }

    TestStatsEqual(const DATA& actual, unsigned count_missing, T missing, T min, double avg, T max)
        : matrix(actual), count_missing(count_missing), missing(missing), min(min), avg(avg), max(max)
    {
    }

    void check(WIBBLE_TEST_LOCPRM) const
    {
        using namespace wibble::tests;
        using namespace std;

        ArrayStats<T> stats;
        bool failed = false;
        if (has_missing)
        {
            stats.fill(missing, matrix);
            if (stats.count_missing != count_missing)
                failed = true;
        } else
            stats.fill(matrix);
        if (stats.min != min) failed = true;
        if (stats.max != max) failed = true;
        if (round(stats.avg * 100) != round(avg*100)) failed = true;

        if (failed)
        {
            std::stringstream ss;
            ss << "stats (";
            if (has_missing)
                ss << "missing: " << stats.count_missing << " ";
            ss << "min: " << (double)stats.min
               << " avg: " << fixed << setprecision(2) << (double)stats.avg
               << " max: " << (double)stats.max
               << ") differ from expected (";
            if (has_missing)
                ss << "missing: " << count_missing << " ";
            ss << "min: " << (double)min
               << " avg: " << fixed << setprecision(2) << (double)avg
               << " max: " << (double)max
               << ")";
            wibble_test_location.fail_test(ss.str());
        }
    }
};

template<typename T>
struct ActualMatrix2D : public wibble::tests::Actual<const cumbac::Matrix2D<T>&>
{
    using wibble::tests::Actual<const cumbac::Matrix2D<T>&>::Actual;

    template<typename... args>
    TestStatsEqual<cumbac::Matrix2D<T>, T> statsEqual(args&&... params) const
    {
        return TestStatsEqual<cumbac::Matrix2D<T>, T>(this->actual, params...);
    }
};

template<typename T>
struct ActualVolume : public wibble::tests::Actual<const cumbac::Volume<T>&>
{
    using wibble::tests::Actual<const cumbac::Volume<T>&>::Actual;

    template<typename... args>
    TestStatsEqual<cumbac::Volume<T>, T> statsEqual(args&&... params) const
    {
        return TestStatsEqual<cumbac::Volume<T>, T>(this->actual, params...);
    }
};

template<typename T>
inline ActualMatrix2D<T> actual(const cumbac::Matrix2D<T>& actual) { return ActualMatrix2D<T>(actual); }

template<typename T>
inline ActualMatrix2D<T> actual(const cumbac::PolarScan<T>& actual) { return ActualMatrix2D<T>(actual); }

template<typename T>
inline ActualMatrix2D<T> actual(const cumbac::Image<T>& actual) { return ActualMatrix2D<T>(actual); }

template<typename T>
inline ActualVolume<T> actual(const cumbac::Volume<T>& actual) { return ActualVolume<T>(actual); }


}

#endif
