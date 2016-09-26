#pragma once

#include<iterator>
#include<algorithm>
#include<iomanip>	// setw()
#include<fstream>

#include"my_matrix.h"

template<typename T>
inline Matrix<T>::Matrix(): rows_(0), colls_(0) {
	std::vector<T>().swap(body_);
}

template<class T>
Matrix<T>::Matrix(unsigned rows, unsigned colls) : rows_(rows), colls_(colls)
{
	// Resize matrix body to contain all elements of matrix
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
inline Matrix<T> Matrix<T>::times_divide(Matrix<T> const & other)
{
	// Check that rows or columns number of right and left matrix are equal
	assert(rows_ == other.rows_, colls_ == other.colls_);

	Matrix<T> result(*this);
#pragma omp parallel for
	for (int i = 0; i < result.body_.size(); ++i)
		result.body_.at(i) /= other.body_.at(i);

	return result;
}

template<typename T>
inline std::pair<unsigned int, unsigned int> Matrix<T>::size() const
{
	return std::make_pair(rows_, colls_);
}

template<typename T>
inline long double Matrix<T>::getSum() const
{
	long double sum{ 0.0 };

#pragma omp parallel for
	for (int i = 0; i < body_.size(); ++i)
		sum += static_cast<long double>(body_.at(i));

	return sum;
}

template<typename T>
inline std::vector<T> Matrix<T>::getRow(unsigned const y) const
{
	// Check that row ID less than number of rows
	assert(y < rows_);
	std::vector<T> result(colls_, T());
	
#pragma omp parallel for
	for (int x = 0; x < colls_; ++x)
		result.at(x) = body_.at(colls_ * y + x);
	return result;
}

template<typename T>
inline void Matrix<T>::setRow(unsigned const y, std::vector<T> const & row)
{
	// Check that std::vector<T> row size is equal to columns number of matrix
	assert(colls_ == row.size());

#pragma omp parallel
	for (int x = 0; x < colls_; ++x)
		body_.at(colls_ * y + x) = row.at(x);
}

template<typename T>
inline std::vector<T> Matrix<T>::getColumn(unsigned const x) const
{
	// Check that coll ID less than number of column
	assert(x < colls_);
	std::vector<T> result(rows_ - 2, T());
	
#pragma omp parallel for
	for (int y = 1; y < rows_ - 1; ++y)
		result.at(y - 1) = body_.at(x + y * colls_);

	return result;
}

template<typename T>
inline void Matrix<T>::setColumn(unsigned const x, std::vector<T> const & coll)
{
	// Check that std::vector<T> coll size is equal to rows number of matrix, bsides 2 (left and right boundary index)
	assert(rows_ == coll.size() + 2);

#pragma omp parallel for
	for (int y = 1; y < rows_ - 1; ++y)
		body_.at(x + y * colls_) = coll.at(y - 1);
}

template<typename T>
inline void Matrix<T>::fillWith(T const value)
{
#pragma omp parallel for
	for (int i = 0; i < body_.size(); ++i)
		body_.at(i) = value;
}

template<typename T>
inline void Matrix<T>::fillWithoughtBoundary(T const value)
{
#pragma omp parallel for
	for (int i = colls_; i < colls_ * (rows_ - 1); ++i) {
		if (i % colls_ == 0 || (i + 1) % colls_ == 0)
			continue;
		body_.at(i) = value;
	}
}

template<typename T>
inline void Matrix<T>::fillColumnWith(int const coll_id, T const value)
{
	// Check that coll ID less than number of column
	assert(coll_id < colls_);

#pragma omp parallel for
	for (int y = 0; y < rows_; ++y)
		body_.at(coll_id + y * colls_) = value;
}

template<typename T>
inline void Matrix<T>::fillRowWith(int const row_id, T const value)
{
	// Check that row ID less than number of rows
	assert(row_id < rows_);

#pragma omp parallel for
	for (int x = 0; x < colls_; ++x)
		body_.at(row_id * colls_ + x) = value;
}

template<typename T>
inline void Matrix<T>::resize(unsigned rows, unsigned colls)
{
	rows_ = rows;
	colls_ = colls;

	// Трюк с освобождением памяти под вектор | body_.clear() Не работает - проверка по body.capasity()
	std::vector<T>().swap(body_);

	// If we swap on matrix size (y,0) or (0, x) we need onlly to allocate memory
	if(rows_ != 0 && colls_ != 0)
		body_.resize(rows_ * colls_, T());
}

template<typename T>
inline void Matrix<T>::writeToFile(std::string value_name, int const time)
{
	using std::endl;

	// Set PATH and NAME of file
	std::string path = "Data/";
	std::string full_name = path + value_name + "[" + std::to_string(rows_) + "x" + std::to_string(colls_) + "]_at_" + std::to_string(time) + "_time_steps.txt";
	
	std::ofstream os;
	os.open(full_name);

	if (!os.good()){
		std::cout << "Error! Can not open file to write.\n";
		return;
	}
	else {
		unsigned cur_position{ 1 };

		for (unsigned y = 0; y < rows_; ++y) {
			for (unsigned x = 0; x < colls_; ++x) {
				os << body_.at(x + y * colls_);
				if (x != colls_ - 1)
					os << ' ';
			}
			os << endl;
		}
	}
	os.close();
}

template<typename T>
inline void Matrix<T>::writeColumnToFile(std::string value_name, int const coll_id, int const time)
{
	// Check that coll_id is less then columns number
	assert(coll_id < colls_);
	using std::endl;

	// Set PATH and NAME of file
	std::string path = "Data/";
	std::string full_name = path + value_name + "_coll(" + std::to_string(coll_id) +")[" + std::to_string(rows_) + "x" + std::to_string(colls_) + "]_at_" + std::to_string(time) + "_time_steps.txt";

	std::ofstream os;
	os.open(full_name);

	if (!os.good()) {
		std::cout << "Error! Can not open file to write.\n";
		return;
	}
	else {
		for (int y = 0; y < rows_; ++y)
			os << body_.at(coll_id + y * colls_) << std::endl;
	}
	os.close();

}

template<typename T>
inline void Matrix<T>::writeRowToFile(std::string value_name, int const row_id, int const time)
{
	// Check that row_id is less then rows number
	assert(row_id < rows_);
	using std::endl;

	// Set PATH and NAME of file
	std::string path = "Data/";
	std::string full_name = path + value_name + "_row(" + std::to_string(row_id) + ")[" + std::to_string(rows_) + "x" + std::to_string(colls_) + "]_at_" + std::to_string(time) + "_time_steps.txt";

	std::ofstream os;
	os.open(full_name);

	if (!os.good()) {
		std::cout << "Error! Can not open file to write.\n";
		return;
	}
	else {
		for (int x = 0; x < colls_; ++x) {
			os << body_.at(row_id * colls_ + x);
			if (x != colls_ - 1)
				os << ' ';
		}
	}
	os.close();
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
