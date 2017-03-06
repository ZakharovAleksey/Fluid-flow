#include "my_matrix_3d.h"
#pragma once

#ifndef MY_MATRIX_3D_IMPL_H
#define MY_MATRIX_3D_IMPL_H

#include<iostream>

template<typename T>
void CompareSize(const Matrix3D<T> & first, const Matrix3D<T> & second)
{
	assert(first.GetDepthNumber() == second.GetDepthNumber() &&
		first.GetRowsNumber() == second.GetRowsNumber() &&
		first.GetCollsNumber() == second.GetCollsNumber()
	);
}




template<typename T>
Matrix3D<T>::Matrix3D() : rows_(0), colls_(0), depth_(0) {}

template<typename T>
inline Matrix3D<T>::Matrix3D(int rows, int colls, int depth) : Matrix2D<T>(rows, colls), depth_(depth)
{
	body_.resize(GetTotalSize(), T());

	for (int i = 0; i < GetTotalSize(); ++i)
		body_.at(i) = rand() % 10;
}


template<typename T>
Matrix3D<T>::~Matrix3D() {}




template<typename T>
inline Matrix3D<T> const Matrix3D<T>::ScalarMultiplication(Matrix3D<T> const & other)
{
	CompareSize(*this, other);

	Matrix3D<T> res(*this);
#pragma omp parallel for
	for (int i = 0; i < GetTotalSize(); ++i)
		res.at(i) *= other.body_.at(i);

	return res;
}

template<typename T>
inline void Matrix3D<T>::Swap(Matrix3D<T>& other)
{
	std::swap(depth_, other.depth_);
	std::swap(rows_, other.rows_);
	std::swap(colls_, other.colls_);

	std::swap(body_, other.body_);
}




template<typename T1>
inline const Matrix3D<T1> operator-(const Matrix3D<T1>& right)
{
	Matrix3D<T1> res(right);
#pragma omp parallel for
	for (int i = 0; i < right.GetTotalSize(); ++i)
		res.body_.at(i) *= -1;

	return res;
}




template<typename T1>
Matrix3D<T1>& operator+=(Matrix3D<T1>& left, const Matrix3D<T1>& right)
{
	CompareSize(left, right);

#pragma omp parallel for
	for (int i = 0; i < left.GetTotalSize(); ++i)
		left.body_.at(i) += right.body_.at(i);

	return left;
}

template<typename T1>
inline Matrix3D<T1>& operator+=(Matrix3D<T1>& left, const T1 & value)
{
#pragma omp parallel for
	for (int i = 0; i < left.GetTotalSize(); ++i)
		left.body_.at(i) += value;

	return left;
}

template<typename T1>
Matrix3D<T1> const operator+(const Matrix3D<T1>& left, const T1 & right)
{
	Matrix3D<T1> res(left);
	return res += right;
}

template<typename T1>
Matrix3D<T1> const operator+(const T1 & left, const Matrix3D<T1>& right)
{
	return right + left;
}

template<typename T>
Matrix3D<T> const operator+(const Matrix3D<T> & left, const Matrix3D<T> & right)
{
	Matrix3D<T> res(left);
	res += right;
	return res;
}




template<typename T1>
inline Matrix3D<T1>& operator-=(Matrix3D<T1>& left, const Matrix3D<T1> & right)
{
	CompareSize(left, right);

#pragma parallel for
	for (int i = 0; i < left.GetTotalSize(); ++i)
		left.body_.at(i) -= right.body_.at(i);

	return left;
}

template<typename T1>
inline Matrix3D<T1>& operator-=(Matrix3D<T1>& left, const T1 & right)
{
#pragma parallel for
	for (int i = 0; i < left.GetTotalSize(); ++i)
		left.body_.at(i) -= right;

	return left;
}

template<typename T1>
inline Matrix3D<T1> const operator-(const Matrix3D<T1>& left, const T1 & right)
{
	Matrix3D<T1> res(left);
	return res -= right;
}

template<typename T1>
inline Matrix3D<T1> const operator-(const T1 & left, const Matrix3D<T1>& right)
{
	return - right + left;
}

template<typename T>
Matrix3D<T> const operator-(const Matrix3D<T> & left, const Matrix3D<T> & right)
{
	Matrix3D<T> res(left);
	res -= right;
	return res;
}

template<typename T1>
inline Matrix3D<T1>& operator*=(Matrix3D<T1>& left, const T1 & right)
{
#pragma omp parallel for
	for (int i = 0; i < left.GetTotalSize(); ++i)
		left.body_.at(i) *= right;

	return left;
}

template<typename T1>
inline const Matrix3D<T1> operator*(const Matrix3D<T1>& left, const T1 & rigth)
{
	Matrix3D<T1> res(left);
	return res *= rigth;
}

template<typename T1>
inline const Matrix3D<T1> operator*(const T1 & left, const Matrix3D<T1>& right)
{
	return right * left;
}




template<typename T1>
inline Matrix3D<T1>& operator/=(Matrix3D<T1>& left, const T1 & right)
{
#pragma omp parallel for
	for (int i = 0; i < left.GetTotalSize(); ++i)
		left.body_.at(i) /= right;

	return left;
}

template<typename T1>
inline const Matrix3D<T1> operator/(const Matrix3D<T1>& left, const T1 & rigth)
{
	Matrix3D<T1> res(left);
	return res /= rigth;
	
}




template<typename T1>
std::ostream & operator<<(std::ostream & os, Matrix3D<T1> const & matrix) {
	using std::endl;

	os.precision(3);

	for (int z = 0; z < matrix.GetDepthNumber(); ++z)
	{
		os << "Depth layer: " << z << endl;
		for (int y = 0; y < matrix.GetRowsNumber(); ++y)
		{
			for (int x = 0; x < matrix.GetCollsNumber(); ++x)
			{
				std::cout << std::setw(7) << matrix(z, y, x);
			}
			os << endl;
		}
		os << endl;
	}

	os << endl;

	return os;
}

#endif // !MY_MATRIX_3D_IMPL_H
