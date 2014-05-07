#ifndef ARCHIVIATORE_TEST_UTILS_CLASS_H
#define ARCHIVIATORE_TEST_UTILS_CLASS_H

#include <volume.h>

namespace testradar {

template<typename T>
struct ArrayStats
{
    bool first = true;
    T min = 0;
    T max = 0;
    double avg = 0;
    unsigned count_zeros = 0;
    unsigned count_ones = 0;

    ArrayStats() {}

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

    void fill(const T* arr, unsigned size)
    {
        for (unsigned i = 0; i < size; ++i)
            this->count_sample(arr[i], size);
    }

    void fill(const cumbac::Matrix2D<T>& arr)
    {
        for (int i = 0; i < arr.rows() * arr.cols(); ++i)
            this->count_sample(arr.data()[i], arr.rows() * arr.cols());
    }

    void fill(const cumbac::PolarScan<T>& arr)
    {
        for (int i = 0; i < arr.rows() * arr.cols(); ++i)
            this->count_sample(arr.data()[i], arr.rows() * arr.cols());
    }

    void fill(const cumbac::Volume<T>& vol)
    {
        for (unsigned i = 0; i < vol.size(); ++i)
            fill(vol.scan(i));
    }

    template<int A, int B>
    void fill2(const T (&arr)[A][B])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
                this->count_sample(arr[i][j], A * B);
    }

    template<int A, int B, int C>
    void fill3(const T (&arr)[A][B][C])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
                for (int k = 0; k < C; ++k)
                    this->count_sample(arr[i][j][k], A * B * C);
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

}

#endif
