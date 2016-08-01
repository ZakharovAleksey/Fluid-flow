#pragma once

#include"distribution_func.h"

template<typename T>
inline DistributionFunction<T>::DistributionFunction() : rows_(0), colls_(0) {
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).resize(rows_, colls_);
}

template<typename T>
DistributionFunction<T>::DistributionFunction(unsigned rows, unsigned colls): rows_(rows), colls_(colls) {

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).resize(rows_, colls_);
}

template<typename T>
DistributionFunction<T>::~DistributionFunction() {}

template<typename T>
inline DistributionFunction<T>::DistributionFunction(DistributionFunction const & other) :
	rows_(other.rows_), colls_(other.colls_)
{
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		dfunc_body_.at(q).resize(rows_, colls_);
		dfunc_body_.at(q) = other.dfunc_body_.at(q);
	}
}

template<typename T>
inline void DistributionFunction<T>::swap(DistributionFunction & dist_func)
{
	std::swap(rows_, dist_func.rows_);
	std::swap(colls_, dist_func.colls_);

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).swap(dist_func.dfunc_body_.at(q));
}

template<typename T>
inline Matrix<T>& DistributionFunction<T>::operator[](unsigned q)
{
	assert(q < kQ);
	return dfunc_body_.at(q);
}

template<typename T>
inline void DistributionFunction<T>::fill(T const value)
{
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).fill_with(value);
}

template<typename T>
inline void DistributionFunction<T>::resize(unsigned rows, unsigned colls)
{
	rows_ = rows;
	colls_ = colls;

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).resize(rows_, colls_);
}

template<typename T>
inline std::array<Matrix<T>, kQ> DistributionFunction<T>::get_values_on_upper_boundary() const
{
	std::array<Matrix<T>, kQ> result;
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		result.at(q).resize(1, colls_);
	}

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		result.at(q) = dfunc_body_.at(q).get_row(1);
		//std::cout << "f[" << q << "] = " << result.at(q);
	}

	return result;
}

template<typename T>
inline std::array<Matrix<T>, kQ> DistributionFunction<T>::get_values_on_bottom_boundary() const
{
	std::array<Matrix<T>, kQ> result;
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		result.at(q).resize(1, colls_);
	}

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		result.at(q) = dfunc_body_.at(q).get_row(rows_ - 2);
		//std::cout << "f[" << q << "] = " << result.at(q);
	}

	return result;
}

template<typename T>
inline std::array<Matrix<T>, kQ> DistributionFunction<T>::get_values_on_left_boundary() const
{
	std::array<Matrix<T>, kQ> result;
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		result.at(q).resize(1, rows_ - 2);
	}

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		result.at(q) = dfunc_body_.at(q).get_coll(1);
		//std::cout << "f[" << q << "] = " << result.at(q);
	}

	return result;
}

template<typename T>
inline std::array<Matrix<T>, kQ> DistributionFunction<T>::get_values_on_right_boundary() const
{
	std::array<Matrix<T>, kQ> result;
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		result.at(q).resize(1, rows_ - 2);
	}

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		result.at(q) = dfunc_body_.at(q).get_coll(colls_ - 2);
		std::cout << "f[" << q << "] = " << result.at(q);
	}

	return result;
}

template<typename T1>
std::ostream & operator<<(std::ostream & os, DistributionFunction<T1> const & dist_func) {
	unsigned i{ 0 };
	for (auto matrix : dist_func.dfunc_body_) {
		os << "Distribution func f[" << i++ << "] ------------------- \n";
		os << matrix;
	}

	return os;
}