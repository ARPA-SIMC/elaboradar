#ifndef ELABORADAR_MATRIX_H
#define ELABORADAR_MATRIX_H

#include <stdexcept>
#include <cmath>
#include <Eigen/Dense>

namespace cumbac {

/// Base for all matrices we use, since we rely on row-major data
template<typename T>
struct Matrix2D : public Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
    using Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Matrix;

    T* row_ptr(unsigned y) { return this->data() + y * this->cols(); }
};

template<typename T>
struct Image : public Matrix2D<T>
{
    Image(unsigned sx, unsigned sy=0)
        : Matrix2D<T>(Matrix2D<T>::Zero(sy ? sy : sx, sx)) {}

};

}

#endif
