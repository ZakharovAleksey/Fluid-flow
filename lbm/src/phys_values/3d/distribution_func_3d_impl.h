#pragma once

#ifndef DISTRIBUTION_FUNC_IMPL_3D_H
#define DISTRIBUTION_FUNC_IMPL_3D_H

#include"distribution_func_3d.h"

template<typename T>
DistributionFunction3D<T>::DistributionFunction3D(int depth, int rows, int colls) : depth_(depth), rows_(rows), colls_(colls) 
{
	for (int q = 0; q < kQ3d; ++q)
	{
		body_.at(q).Resize(rows_, colls_, depth_);
	}
}

template<typename T>
DistributionFunction3D<T>::~DistributionFunction3D() {}



template<typename T>
inline Matrix3D<T>& DistributionFunction3D<T>::operator[](const int q)
{
	return body_.at(q);
}

template<typename T>
inline const Matrix3D<T>& DistributionFunction3D<T>::operator[](const int q) const
{
	return body_.at(q);
}

template<typename T>
inline void DistributionFunction3D<T>::Swap(DistributionFunction3D<T>& other)
{
	std::swap(depth_, other.depth_);
	std::swap(rows_, other.rows_);
	std::swap(colls_, other.colls_);

	std::swap(body_, other.body_);
}

template<typename T>
inline std::vector<T> DistributionFunction3D<T>::GetTopBoundaryValues(int const q) const
{
	return body_.at(q).GetRow(1);
}

template<typename T>
inline std::vector<T> DistributionFunction3D<T>::GetBottomBoundaryValue(int const q) const
{
	return body_.at(q).GetRow(rows_ - 2);
}

template<typename T>
inline std::vector<T> DistributionFunction3D<T>::GetLeftBoundaryValue(int const q) const
{
	return body_.at(q).GetColumn(1);
}

template<typename T>
inline std::vector<T> DistributionFunction3D<T>::GetRightBoundaryValue(int const q) const
{
	return body_.at(q).GetColumn(colls_ - 2);
}

template<typename T>
inline void DistributionFunction3D<T>::SetTopBoundaryValue(int const q, std::vector<T> const & row)
{
	body_.at(q).SetRow(1, row);
}

template<typename T>
inline void DistributionFunction3D<T>::SetBottomBoundaryValue(int const q, std::vector<T> const & row)
{
	body_.at(q).SetRow(rows_ - 2, row);
}

template<typename T>
inline void DistributionFunction3D<T>::SetLeftBoundaryValue(int const q, std::vector<T> const & coll)
{
	body_.at(q).SetColumn(1, coll);
}

template<typename T>
inline void DistributionFunction3D<T>::SetRightBoundaryValue(int const q, std::vector<T> const & coll)
{
	body_.at(q).SetColumn(colls_ - 2, coll);
}


template<typename T1>
inline const DistributionFunction3D<T1> operator-(const DistributionFunction3D<T1>& right)
{
	DistributionFunction3D<T1> res(right);
	for (int q = 0; q < kQ3d; ++q)
		res[q] *= -1;

	return res;
}

template<typename T1>
inline DistributionFunction3D<T1>& operator+=(DistributionFunction3D<T1>& left, const DistributionFunction3D<T1>& right)
{
	for (int q = 0; q < kQ3d; ++q)
		left[q] += right[q];

	return left;
}

template<typename T1>
inline DistributionFunction3D<T1>& operator+=(DistributionFunction3D<T1>& left, const T1 & right)
{
	for (int q = 0; q < kQ3d; ++q)
		left[q] += right;

	return left;
}

template<typename T1>
inline DistributionFunction3D<T1>& operator-=(DistributionFunction3D<T1>& left, const DistributionFunction3D<T1>& right)
{
	for (int q = 0; q < kQ3d; ++q)
		left[q] -= right[q];

	return left;
}

template<typename T1>
inline DistributionFunction3D<T1>& operator-=(DistributionFunction3D<T1>& left, const T1 & right)
{
	for (int q = 0; q < kQ3d; ++q)
		left[q] -= right;

	return left;
}

template<typename T1>
inline DistributionFunction3D<T1>& operator*=(DistributionFunction3D<T1>& left, const T1 & right)
{
	for (int q = 0; q < kQ3d; ++q)
		left[q] *= right;

	return left;
}

template<typename T1>
inline DistributionFunction3D<T1>& operator/=(DistributionFunction3D<T1>& left, const T1 & right)
{
	for (int q = 0; q < kQ3d; ++q)
		left[q] /= right;

	return left;
}




template<typename T>
const DistributionFunction3D<T> operator+(const DistributionFunction3D<T>& left, const DistributionFunction3D<T>& right)
{
	DistributionFunction3D<T> res(left);
	return res += right;
}

template<typename T>
const DistributionFunction3D<T> operator+(const DistributionFunction3D<T>& left, const T & right)
{
	DistributionFunction3D<T> res(left);
	return res += right;
}

template<typename T>
const DistributionFunction3D<T> operator+(const T & left, const DistributionFunction3D<T>& right)
{
	return right + left;
}


template<typename T>
const DistributionFunction3D<T> operator-(const DistributionFunction3D<T>& left, const DistributionFunction3D<T>& right)
{
	DistributionFunction3D<T> res(left);
	return res -= right;
}

template<typename T>
const DistributionFunction3D<T> operator-(const DistributionFunction3D<T>& left, const T & right)
{
	DistributionFunction3D<T> res(left);
	return res -= right;
}

template<typename T>
const DistributionFunction3D<T> operator-(const T & left, const DistributionFunction3D<T>& right)
{
	DistributionFunction3D<T> res(right);
	return - right + left;
}

template<typename T>
const DistributionFunction3D<T> operator*(const DistributionFunction3D<T>& left, const T & right)
{
	DistributionFunction3D<T> res(left);
	return res *= right;
}


template<typename T>
const DistributionFunction3D<T> operator*(const T & left, const DistributionFunction3D<T>& right)
{
	return right * left;
}

template<typename T>
const DistributionFunction3D<T> operator/(const DistributionFunction3D<T>& left, const T & right)
{
	DistributionFunction3D<T> res(left);
	return res /= right;
}


template<typename T1>
inline std::ostream & operator<<(std::ostream & os, const DistributionFunction3D<T1>& dist_func)
{
	using std::endl;
	os.precision(3);

	for (int q = 0; q < kQ3d; ++q)
	{
		os << q << " component of 3d distribution function:" << endl;
		os << dist_func.body_.at(q);
	}

	return os;
}

#endif // !DISTRIBUTION_FUNC_IMPL_3D_H
