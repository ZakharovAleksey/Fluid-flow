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
inline void DistributionFunction<T>::fill_boundaries(T const value)
{
//#pragma omp parallel for
	for (int q = 1; q < kQ; ++q) {	// Начинаем с 1 так как f[0] никуда не двигается
		switch (q)
		{
		case 1:
			dfunc_body_.at(1).fill_coll_with(colls_ - 1, value);
			break;
		case 2:
			dfunc_body_.at(2).fill_row_with(0, value);
			break;
		case 3:
			dfunc_body_.at(3).fill_coll_with(0, value);
			break;
		case 4:
			dfunc_body_.at(4).fill_row_with(rows_ - 1, value);
			break;
		case 5:
			dfunc_body_.at(5).fill_coll_with(colls_ - 1, value);
			dfunc_body_.at(5).fill_row_with(0, value);
			break;
		case 6:
			dfunc_body_.at(6).fill_coll_with(0, value);
			dfunc_body_.at(6).fill_row_with(0, value);
			break;
		case 7:
			dfunc_body_.at(7).fill_coll_with(0, value);
			dfunc_body_.at(7).fill_row_with(rows_ - 1, value);
			break;
		case 8:
			dfunc_body_.at(8).fill_coll_with(colls_ - 1, value);
			dfunc_body_.at(8).fill_row_with(rows_ - 1, value);
			break;
		default:
			break;
		}
	}
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
inline std::pair<unsigned int, unsigned int> DistributionFunction<T>::size() const
{
	return std::make_pair(rows_, colls_);
}

template<typename T>
inline std::vector<T> DistributionFunction<T>::get_top_boundary_val(int const q) const
{
	return dfunc_body_.at(q).get_row(1);
}

template<typename T>
inline std::vector<T> DistributionFunction<T>::get_bottom_boundary_val(int const q) const
{
	return dfunc_body_.at(q).get_row(rows_ - 2);
}

template<typename T>
inline std::vector<T> DistributionFunction<T>::get_left_boundary_val(int const q) const
{
	return dfunc_body_.at(q).get_coll(1);
}

template<typename T>
inline std::vector<T> DistributionFunction<T>::get_right_boundary_val(int const q) const
{
	return dfunc_body_.at(q).get_coll(colls_ - 2);
}

template<typename T>
inline void DistributionFunction<T>::set_top_boundary_value(int const q, std::vector<T> const & row)
{
	dfunc_body_.at(q).set_row(1, row);
}

template<typename T>
inline void DistributionFunction<T>::set_bottom_boundary_value(int const q, std::vector<T> const & row)
{
	dfunc_body_.at(q).set_row(rows_ - 2, row);
}

template<typename T>
inline void DistributionFunction<T>::set_left_boundary_value(int const q, std::vector<T> const & coll)
{
	dfunc_body_.at(q).set_coll(1, coll);
}

template<typename T>
inline void DistributionFunction<T>::set_right_boundary_value(int const q, std::vector<T> const & coll)
{
	dfunc_body_.at(q).set_coll(colls_ - 2, coll);
}

template<typename T>
inline MacroscopicParam<T> DistributionFunction<T>::get_density() const
{
	MacroscopicParam<T> result(rows_, colls_);
	result.fill_with(0.0);

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q)
		result += dfunc_body_.at(q);

	return result;
}

template<typename T>
inline MacroscopicParam<T> DistributionFunction<T>::get_velocity(const double mas[kQ], MacroscopicParam<T> const & density) const
{
	MacroscopicParam<T> result(rows_, colls_);
	result.fill_with(0.0);

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q)
		result += dfunc_body_.at(q) * mas[q];

	// Переписать сразу в return, тк лишнее копирование НО выпадает error
	result.times_divide(density);
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