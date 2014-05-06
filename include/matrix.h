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

    T* row(unsigned y) { return this->data() + y * this->cols(); }
};

}

#endif
