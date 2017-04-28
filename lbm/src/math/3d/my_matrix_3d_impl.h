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
Matrix3D<T>::Matrix3D() : depth_(0), rows_(0), colls_(0) {}

template<typename T>
inline Matrix3D<T>::Matrix3D(int depth, int rows, int colls) : depth_(depth), rows_(rows), colls_(colls)
{
	body_.resize(GetTotalSize(), T());

//	for (int i = 0; i < GetTotalSize(); ++i)
//		body_.at(i) = rand() % 10;
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
		res.body_.at(i) *= other.body_.at(i);

	return res;
}

template<typename T>
inline Matrix3D<T> const Matrix3D<T>::TimesDivide(Matrix3D<T> const & other)
{
	CompareSize(*this, other);

	Matrix3D<T> res(*this);
#pragma omp parallel for
	for (int i = 0; i < GetTotalSize(); ++i)
		res.body_.at(i) /= other.body_.at(i);

	return res;
}

template<typename T>
inline void Matrix3D<T>::Resize(int new_rows_numb, int new_colls_numb, int new_depth_numb = 0)
{
	rows_ = new_rows_numb;
	colls_ = new_colls_numb;
	depth_ = new_depth_numb;

	// Трюк с освобождением памяти под вектор | body_.clear() Не работает - проверка по body.capasity()
	std::vector<T>().swap(body_);

	// If we swap on matrix size (y,0) or (0, x) we need onlly to allocate memory
	if (depth_ != 0 && rows_ != 0 && colls_ != 0)
		body_.resize(depth_ * rows_ * colls_, T());
}

template<typename T>
inline long double Matrix3D<T>::GetSum() const
{
	long double sum{ 0.0 };

#pragma omp parallel for
	for (int i = 0; i < GetTotalSize(); ++i)
	{
		#pragma omp atomic
		sum += static_cast<long double>(body_.at(i));
	}

	return sum;


}

template<typename T>
inline std::vector<T> Matrix3D<T>::GetRow(unsigned const y) const
{
	// Check that row ID less than number of rows
	assert(y < rows_);
	std::vector<T> result(depth_ * colls_, T());

	for (int z = 0; z < depth_; ++z)
	{
		#pragma omp parallel for
		for (int x = 0; x < colls_; ++x)
			result.at(z * colls_ + x) = this->operator()(z, y, x);
	}
	return result;
}

template<typename T>
inline void Matrix3D<T>::SetRow(unsigned const y, std::vector<T> const & row)
{
	// Check that std::vector<T> row size is equal to columns number of matrix
	assert(colls_ * depth_ == row.size());

	for (int z = 0; z < depth_; ++z)
	{
		#pragma omp parallel
		for (int x = 0; x < colls_; ++x)
			this->operator()(z,y,x) = row.at(z * colls_ + x);
	}
}

template<typename T>
inline std::vector<T> Matrix3D<T>::GetColumn(unsigned const x) const
{
	// Check that coll ID less than number of column
	assert(x < colls_);
	std::vector<T> result(depth_ * (rows_ - 2), T());

	int curId = 0;
	for (int z = 0; z < depth_; ++z)
	{
		#pragma omp parallel for
		for (int y = 1; y < rows_ - 1; ++y)
			result.at(curId++) = this->operator()(z, y, x);
			//result.at(y - 1) = body_.at(x + y * colls_);
	}
	return result;
}

template<typename T>
inline void Matrix3D<T>::SetColumn(unsigned const x, std::vector<T> const & coll)
{
	// Check that std::vector<T> coll size is equal to rows number of matrix, bsides 2 (left and right boundary index)
	assert(depth_ * (rows_ - 2) == coll.size());

	int curId = 0;
	for (int z = 0; z < depth_; ++z)
	{
		#pragma omp parallel for
		for (int y = 1; y < rows_ - 1; ++y)
			this->operator()(z, y, x) = coll.at(curId++);
	}
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
