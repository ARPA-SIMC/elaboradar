#ifndef ELABORADAR_MATRIX_H
#define ELABORADAR_MATRIX_H

#include <stdexcept>
#include <cmath>
#include <Eigen/Dense>

namespace cumbac {

/// Generic 2D matrix
template<typename T>
struct Matrix2D : Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
    // Create a new matrix, with all elements set to a fill value
    Matrix2D(unsigned SY, unsigned int SX, const T& fill=0)
        : Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(
                Matrix2D::Constant(SY, SX, fill)) {}

    // Number of elements in the matrix
    size_t size() const { return this->rows() * this->cols(); }

    T* row(unsigned y) { return this->data() + y * this->cols(); }

    T min() const { return this->minCoeff(); }
    T max() const { return this->maxCoeff(); }
    T avg() const
    {
        unsigned size = this->rows() * this->cols();
        double mean;
        for (unsigned i = 0; i < size; ++i)
            mean += (double)this->data()[i] / size;
        return round(mean);
    }
};

}

#endif
