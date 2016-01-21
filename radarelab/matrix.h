/**
 *  @file
 *  @ingroup radarelab 
*/
#ifndef RADARELAB_MATRIX_H
#define RADARELAB_MATRIX_H

#include <stdexcept>
#include <cmath>
#include <Eigen/Dense>

#include <iostream>

namespace radarelab {

template<typename Scalar> struct Exponential {
EIGEN_EMPTY_STRUCT_CTOR(Exponential)
//typedef Scalar result_type;
Scalar operator()(const Scalar& base, const Scalar& exponent) const { return std::pow(base,exponent); }
};

template<typename Scalar> struct Logarithm {
EIGEN_EMPTY_STRUCT_CTOR(Logarithm)
//typedef Scalar result_type;
Scalar operator()(const Scalar& base, const Scalar& argument) const { return std::log(argument)/std::log(base); }
};

template<typename Scalar> struct Log10 {
EIGEN_EMPTY_STRUCT_CTOR(Log10)
//typedef Scalar result_type;
Scalar operator()(const Scalar& argument) const { return std::log10(argument); }
};

/// Base for all matrices we use, since we rely on row-major data
template<class T>
struct Matrix2D : public Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
    using Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Matrix;

    T* row_ptr(unsigned y) { return this->data() + y * this->cols(); }

    Matrix2D<T> log10()
    {
        return this->unaryExpr(Log10<T>());
    }
    Matrix2D<T> exp10()
    {
        Matrix2D<T> base(Matrix2D<T>::Constant(this->rows(),this->cols(),10.));
        return base.binaryExpr(*this,Exponential<T>());
    }


};

template<typename T>
struct Image : public Matrix2D<T>
{
    Image(unsigned sx, unsigned sy=0)
        : Matrix2D<T>(Matrix2D<T>::Zero(sy ? sy : sx, sx)) {}

};

}

#endif
