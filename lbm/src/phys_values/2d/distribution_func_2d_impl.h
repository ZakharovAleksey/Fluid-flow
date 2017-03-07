#pragma once

#include"distribution_func_2d.h"

template<typename T>
inline DistributionFunction<T>::DistributionFunction() : rows_(0), colls_(0) 
{
	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).Resize(rows_, colls_);
}

template<typename T>
DistributionFunction<T>::DistributionFunction(unsigned rows, unsigned colls): rows_(rows), colls_(colls) 
{
	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).resize(rows_, colls_);
}

template<typename T>
DistributionFunction<T>::~DistributionFunction() {}

template<typename T>
inline DistributionFunction<T>::DistributionFunction(DistributionFunction const & other) :
	rows_(other.rows_), colls_(other.colls_)
{
	for (int q = 0; q < kQ; ++q) 
	{
		dfunc_body_.at(q).resize(rows_, colls_);
		dfunc_body_.at(q) = other.dfunc_body_.at(q);
	}
}

template<typename T>
inline void DistributionFunction<T>::swap(DistributionFunction & dist_func)
{
	std::swap(rows_, dist_func.rows_);
	std::swap(colls_, dist_func.colls_);

	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).swap(dist_func.dfunc_body_.at(q));
}

template<typename T>
inline Matrix2D<T>& DistributionFunction<T>::operator[](unsigned q)
{
	assert(q < kQ);
	return dfunc_body_.at(q);
}

template<typename T>
inline void DistributionFunction<T>::fillWithoutBoundaries(T const value)
{
	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).FillWith(value);
}

template<typename T>
inline void DistributionFunction<T>::fillBoundaries(T const value)
{
	for (int q = 1; q < kQ; ++q) {	// Начинаем с 1 так как f[0] никуда не двигается
		switch (q)
		{
		case 1:
			dfunc_body_.at(1).FillColumnWith(colls_ - 1, value);
			break;
		case 2:
			dfunc_body_.at(2).FillRowWith(0, value);
			break;
		case 3:
			dfunc_body_.at(3).FillColumnWith(0, value);
			break;
		case 4:
			dfunc_body_.at(4).FillRowWith(rows_ - 1, value);
			break;
		case 5:
			dfunc_body_.at(5).FillColumnWith(colls_ - 1, value);
			dfunc_body_.at(5).FillRowWith(0, value);
			break;
		case 6:
			dfunc_body_.at(6).FillColumnWith(0, value);
			dfunc_body_.at(6).FillRowWith(0, value);
			break;
		case 7:
			dfunc_body_.at(7).FillColumnWith(0, value);
			dfunc_body_.at(7).FillRowWith(rows_ - 1, value);
			break;
		case 8:
			dfunc_body_.at(8).FillColumnWith(colls_ - 1, value);
			dfunc_body_.at(8).FillRowWith(rows_ - 1, value);
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

	for (int q = 0; q < kQ; ++q)
		dfunc_body_.at(q).Resize(rows_, colls_);
}

template<typename T>
inline std::pair<unsigned int, unsigned int> DistributionFunction<T>::size() const
{
	return std::make_pair(rows_, colls_);
}

template<typename T>
inline std::vector<T> DistributionFunction<T>::getTopBoundaryValues(int const q) const
{
	return dfunc_body_.at(q).GetRow(1);
}

template<typename T>
inline std::vector<T> DistributionFunction<T>::getBottomBoundaryValue(int const q) const
{
	return dfunc_body_.at(q).GetRow(rows_ - 2);
}

template<typename T>
inline std::vector<T> DistributionFunction<T>::getLeftBoundaryValue(int const q) const
{
	return dfunc_body_.at(q).GetColumn(1);
}

template<typename T>
inline std::vector<T> DistributionFunction<T>::getRightBoundaryValue(int const q) const
{
	return dfunc_body_.at(q).GetColumn(colls_ - 2);
}

template<typename T>
inline void DistributionFunction<T>::setTopBoundaryValue(int const q, std::vector<T> const & row)
{
	dfunc_body_.at(q).SetRow(1, row);
}

template<typename T>
inline void DistributionFunction<T>::setBottomBoundaryValue(int const q, std::vector<T> const & row)
{
	dfunc_body_.at(q).SetRow(rows_ - 2, row);
}

template<typename T>
inline void DistributionFunction<T>::setLeftBoundaryValue(int const q, std::vector<T> const & coll)
{
	dfunc_body_.at(q).SetColumn(1, coll);
}

template<typename T>
inline void DistributionFunction<T>::setRightBoundaryValue(int const q, std::vector<T> const & coll)
{
	dfunc_body_.at(q).SetColumn(colls_ - 2, coll);
}

template<typename T>
inline MacroscopicParam<T> DistributionFunction<T>::calculateDensity() const
{
	MacroscopicParam<T> result(rows_, colls_);
	result.FillWith(0.0);

	for (int q = 0; q < kQ; ++q)
		result += dfunc_body_.at(q);

	return result;
}

template<typename T>
inline MacroscopicParam<T> DistributionFunction<T>::calculateVelocity(const double mas[kQ], MacroscopicParam<T> const & density) const
{
	MacroscopicParam<T> result(rows_, colls_);
	result.FillWith(0.0);

	for (int q = 0; q < kQ; ++q)
		result += dfunc_body_.at(q) * mas[q];

	// Переписать сразу в return, тк лишнее копирование НО выпадает error
	result.TimesDivide(density);
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