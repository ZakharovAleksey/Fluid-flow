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
inline void Matrix3D<T>::FillBoundarySideWalls(const T value)
{
	for (int z = 0; z < depth_; ++z)
	{
		for (int y = 0; y < rows_; ++y)
		{
			this->operator()(z, y, 0) = value;
			this->operator()(z, y, colls_ - 1) = value;
		}

		for (int x = 0; x < rows_; ++x)
		{
			this->operator()(z, 0, x) = value;
			this->operator()(z, rows_ - 1, x) = value;
		}
	}
}

template<typename T>
inline void Matrix3D<T>::FillLayer(const int z, const T value)
{
	for (int y = 0; y < rows_; ++y)
		for (int x = 0; x < colls_; ++x)
		{
			this->operator()(z, y, x) = value;
		}
}

template<typename T>
inline void Matrix3D<T>::FillTopBottomWalls(const T value)
{
	FillLayer(0, value);
	FillLayer(depth_ - 1, value);
}

template<typename T>
inline void Matrix3D<T>::FillWithoutBoundary(const T value)
{
	for (int z = 1; z < depth_ - 1; ++z)
		for (int y = 1; y < rows_ - 1; ++y)
			for (int x = 1; x < colls_ - 1; ++x)
				body_.at(z * rows_ * colls_ + y * colls_ + x) = value;

}

template<typename T>
inline void Matrix3D<T>::FillWith(const T value)
{
	for (auto & i : body_)
		i = value;
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
inline std::vector<T> Matrix3D<T>::GetTBLayer(const int z) const
{
	// We could take layers with are NOT boundaries
	assert(z > 0 && z < depth_ - 1);
	// This layer is FULL : 
	// all Oxy plane including all elements
	std::vector<T> res(rows_ * colls_, T());

	for (int y = 0; y < rows_; ++y)
	{
		for (int x = 0; x < colls_; ++x)
			res.at(y * colls_ + x) = this->operator()(z, y, x);
	}

	return res;
}

template<typename T>
inline std::vector<T> Matrix3D<T>::GetLRLayer(const int x) const
{
	// We could take layers with are NOT boundaries
	assert(x > 0 && x < colls_ - 1);
	// This layer is NOT FULL : upper and bottom elements is belongs to TOP and BOTTOM BS respectively
	// all elements of Oyz expect upper and lower elements
	std::vector<T> res(rows_ * (depth_ - 2), T());

	for (int z = 1; z < depth_ - 1; ++z)
	{
		for(int y = 0; y < rows_; ++y)
			res.at((z - 1) * (rows_) + y) = this->operator()(z, y, x); // (z-1) to fit in vector range
	}

	return res;
}

template<typename T>
inline std::vector<T> Matrix3D<T>::GetNFLayer(const int y) const
{
	// We could take layers with are NOT boundaries
	assert(y > 0 && y < rows_ - 1);
	// This layer is NOT FULL : upper, bottom, left and right elements is belongs to T B R L boundaries respectively
	// all elements of Ozx expect upper, lower, left and right elements
	std::vector<T> res((colls_ - 2) * (depth_ - 2), T());

	for (int z = 1; z < depth_ - 1; ++z)
	{
		for (int x = 1; x < colls_ - 1; ++x)
			res.at((z - 1) * (colls_ - 2) + (x - 1)) = this->operator()(z, y, x); // (z-1),(x-1),(colls-2) to fit in vector range
	}

	return res;
}

template<typename T>
inline void Matrix3D<T>::SetTBLayer(unsigned const z, std::vector<T> const & layer)
{
	assert(layer.size() == rows_ * colls_);

	for (int y = 0; y < rows_; ++y)
	{
		for (int x = 0; x < colls_; ++x)
			this->operator()(z, y, x) = layer.at(y * colls_ + x);
	}
	
}

template<typename T>
inline void Matrix3D<T>::SetLRLayer(unsigned const x, std::vector<T> const & layer)
{
	assert(layer.size() == rows_ * (depth_ - 2));
	
	for (int z = 1; z < depth_ - 1; ++z)
	{
		for (int y = 0; y < rows_; ++y)
			this->operator()(z, y, x) = layer.at((z - 1) * (rows_) + y);
	}
}

template<typename T>
inline void Matrix3D<T>::SetNFLayer(unsigned const y, std::vector<T> const & layer)
{
	assert(layer.size() == (colls_ - 2) * (depth_ - 2));

	for (int z = 1; z < depth_ - 1; ++z)
	{
		for (int x = 1; x < colls_ - 1; ++x)
			this->operator()(z, y, x) = layer.at((z - 1) * (colls_ - 2) + (x - 1));
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
inline const Matrix3D<T1> operator-(const Matrix3D<T1>& right)
{
	Matrix3D<T1> res(right);
#pragma omp parallel for
	for (int i = 0; i < right.GetTotalSize(); ++i)
		res.body_.at(i) *= -1;

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
