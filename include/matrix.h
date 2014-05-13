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

}

#endif
