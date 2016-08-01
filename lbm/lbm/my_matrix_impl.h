#pragma once

#include<iterator>
#include<algorithm>
#include<iomanip>	// setw()

#include"my_matrix.h"

template<typename T>
inline Matrix<T>::Matrix(): rows_(0), colls_(0) {
	std::vector<T>().swap(body_);
}

template<class T>
Matrix<T>::Matrix(unsigned rows, unsigned colls) : rows_(rows), colls_(colls)
{
	// Если заменяем на матрицу размера (y,0) или (0, x) то освобождаем память.
	if (rows == 0 || colls == 0)
		std::vector<T>().swap(body_);
	else
		body_.resize(rows_ * colls_, 0.0);

	#ifndef TEST_INCLUDE
	for (auto & i : body_)
		i = rand() % 10;
	#endif // !TEST_INCLUDE
}

template<class T>
Matrix<T>::~Matrix() {}

template<typename T>
inline Matrix<T>::Matrix(Matrix<T> const & other) : rows_(other.rows_), colls_(other.colls_)
{
	body_.resize(rows_ * colls_);
	body_ = other.body_;
}

template<typename T>
inline void Matrix<T>::swap(Matrix<T> & other)
{
	using std::swap;
	swap(rows_, other.rows_);
	swap(colls_, other.colls_);

	other.body_.swap(body_);
}


template<typename T1>
inline Matrix<T1> operator+(T1 const left, Matrix<T1> const & right)
{
	return Matrix<T1>(right + left);
}

template<typename T1>
inline Matrix<T1> operator*(T1 const left, Matrix<T1> const & right)
{
	return Matrix<T1>(right * left);
}

template<typename T>
inline Matrix<T> Matrix<T>::scalar_multiplication(Matrix<T> const & other)
{
	Matrix<T> result(*this);
#pragma omp parallel for
	for (int i = 0; i < rows_ * colls_; ++i)
		result.body_.at(i) *= other.body_.at(i);

	return result;
}

template<typename T>
inline std::pair<unsigned int, unsigned int> Matrix<T>::size() const
{
	return std::make_pair(rows_, colls_);
}

template<typename T>
inline long double Matrix<T>::get_sum() const
{
	long double sum{ 0.0 };
#pragma omp parallel for
	for (int i = 0; i < body_.size(); ++i)
		sum += static_cast<long double>(body_.at(i));

	return sum;
}

template<typename T>
inline Matrix<T> Matrix<T>::get_row(unsigned const y) const
{
	assert(y < rows_);
	Matrix<T> result(1, colls_);
#pragma omp parallel for
	for (int x = 0; x < colls_; ++x)
		result(0, x) = body_.at(colls_ * y + x);
	return result;
}

template<typename T>
inline void Matrix<T>::set_row(unsigned const y, Matrix<T> const & row)
{
	assert(colls_ == row.size().second && y < rows_ && row.size().first == 1);
#pragma omp parallel
	for (int x = 0; x < colls_; ++x)
		body_.at(colls_ * y + x) = row(0, x);
}

template<typename T>
inline Matrix<T> Matrix<T>::get_coll(unsigned const x) const
{
	assert(x < colls_);
	Matrix<T> result(1, rows_ - 2);
#pragma omp parallel for
	for (int y = 1; y < rows_ - 1; ++y)
		result(0, y - 1) = body_.at(x + y * colls_);

	return result;
}

template<typename T>
inline void Matrix<T>::set_coll(unsigned const x, Matrix<T> const & coll)
{
	assert(rows_ == coll.size().second + 2 && x < colls_ && coll.size().first == 1);
#pragma omp parallel for
	for (int y = 1; y < rows_ - 1; ++y)
		body_.at(x + y * colls_) = coll(0, y - 1);
}

template<typename T>
inline void Matrix<T>::fill_with(T const value)
{
#pragma omp parallel for
	for (int i = 0; i < body_.size(); ++i)
		body_.at(i) = value;
}

template<typename T>
inline void Matrix<T>::fill_withought_boundary(T const value)
{
#pragma omp parallel for
	for (int i = colls_; i < colls_ * (rows_ - 1); ++i) {
		if (i % colls_ == 0 || (i + 1) % colls_ == 0)
			continue;
		body_.at(i) = value;
	}
}

template<typename T>
inline void Matrix<T>::resize(unsigned rows, unsigned colls)
{
	rows_ = rows;
	colls_ = colls;

	// Трюк с освобождением памяти под вектор | body_.clear() Не работает - проверка по body.capasity()
	std::vector<T>().swap(body_);

	body_.resize(rows_ * colls_, T());
}

template<typename T1>
std::ostream & operator<<(std::ostream & os, Matrix<T1> const & matrix) {
	using std::endl;

	os.precision(3);
	unsigned cur_position{ 1 };

	for (unsigned y = 0; y != matrix.rows_; ++y){
		for (unsigned x = 0; x != matrix.colls_; ++x)
			os << std::setw(7) << matrix.body_[x + y * matrix.colls_];
		os << endl;
	}
	


	//std::copy(matrix.body_.begin(), matrix.body_.end(), std::ostream_iterator<T1>(os, " "));
	os << endl;
	
	return os;
}
