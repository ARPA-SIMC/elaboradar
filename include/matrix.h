#ifndef ELABORADAR_MATRIX_H
#define ELABORADAR_MATRIX_H

#include <stdexcept>
#include <cmath>

namespace cumbac {

/// Generic 2D matrix
template<typename T>
struct Matrix2D
{
    const unsigned SX;
    const unsigned SY;
    T* data;

    // Copy constructor, making a copy of the whole matrix
    Matrix2D(const Matrix2D& m)
        : SX(m.SX), SY(m.SY), data(new T[SX * SY])

    {
        for (unsigned i = 0; i < SX * SY; ++i)
            data[i] = m.data[i];
    }

    // Create a new matrix, with all elements set to a fill value
    Matrix2D(unsigned SX, unsigned int SY, const T& fill=0)
        : SX(SX), SY(SY), data(new T[SX * SY])
    {
        for (unsigned i = 0; i < SX * SY; ++i)
            data[i] = fill;
    }

    ~Matrix2D()
    {
        delete[] data;
    }

    // Number of elements in the matrix
    size_t size() const { return SX * SY; }

    // Access a matrix row
    T* operator[](unsigned y) { return data + y * SX; }
    const T* operator[](unsigned y) const { return data + y * SX; }

    Matrix2D& operator=(const Matrix2D& m)
    {
        if (SX != m.SX or SY != m.SY)
            throw std::runtime_error("Matrix2D size mismatch");
        for (unsigned i = 0; i < SX * SY; ++i)
            data[i] = m.data[i];
        return *this;
    }

    T min() const
    {
        T res = data[0];
        for (unsigned i = 0; i < SX * SY; ++i)
            if (data[i] < res)
                res = data[i];
        return res;
    }

    T max() const
    {
        T res = data[0];
        for (unsigned i = 0; i < SX * SY; ++i)
            if (data[i] > res)
                res = data[i];
        return res;
    }

    T avg() const
    {
        double res = 0;
        for (unsigned i = 0; i < SX * SY; ++i)
            res += (double)data[i] / (double)(SX * SY);
        return (T)round(res);
    }
};

}

#endif