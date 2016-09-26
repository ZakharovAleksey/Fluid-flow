#pragma once

#include<iostream>
#include<ostream>
#include <cassert> // setw()
#include<vector>
#include <functional>
#include<string>
#include<omp.h>

#define TEST_INCLUDE 

#ifndef TEST_INCLUDE 
 	#include<cstdlib>
 	#include<ctime>
#endif // !TEST_INCLUDE


/*
	Basic class to this project.  Represent simple matrix object.

	All physical values (density, velocity, distribution function) will be implemented on the basis of this class.
	Overloading implemetation only for operations witch we need in LBM implementation.
*/
template<typename T>
class Matrix
{
public:

#pragma region Constructor

	Matrix();
	Matrix(unsigned rows, unsigned colls);
	virtual ~Matrix();

	Matrix(Matrix<T> const & other);	// Проверить все ли хорошо при вызове функции отсюда

	// Swap left and right matrixes
	void swap(Matrix<T> & other);

#pragma endregion

	Matrix<T> & operator=(Matrix<T> const & other) 
	{
		// Check that rows or colls number of right and left matrix are equal
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		if (this != & other) 
		{
			Matrix<T> temp(other);
			temp.swap(*this);
		}

		return *this;
	}

#pragma region Operator +, += and its overloading

	Matrix<T> & operator+=(Matrix<T> const & other) {
		assert(rows_ == other.rows_ && colls_ == other.colls_);
	#pragma omp parallel for
			for (int i = 0; i < body_.size(); ++i)
				body_.at(i) += other.body_.at(i);

		return *this;
	}

	//! Add to each left matrix element value other
	Matrix<T> & operator+=(T const other) {
		std::for_each(body_.begin(), body_.end(),
			[&](T & value) {value += other; });

		return *this;
	}

	//! Add right matrix to left matrix
	Matrix<T> operator+(Matrix<T> const & other) const {
		// Check that rows or columns number of right and left matrix are equal
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		Matrix<T> result(*this);
	#pragma omp parallel for
		for (int i = 0; i < result.body_.size(); ++i)
			result.body_.at(i) += other.body_.at(i);

		return result;
	}

	//! Add to each left matrix element value other
	Matrix<T> operator+(T const other) const {
		Matrix<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value += other; });

		return result;
	}

	//! Add to each element of right matrix element with value left
	template<typename T1>
	friend Matrix<T1> operator+(T1 const left, Matrix<T1> const & right);

#pragma endregion

#pragma region Operator -, -= and its overloading

	//! Subtraction right matrix from left matrix
	Matrix<T> & operator-=(Matrix<T> const & other) {
		// Check that rows or columns number of right and left matrix are equal
		assert(rows_ == other.rows_ && colls_ == other.colls_);
	#pragma omp parallel for
			for (int i = 0; i < body_.size(); ++i)
				body_.at(i) -= other.body_.at(i);

		return *this;
	}

	//! Subtraction right value from each element of left matrix
	Matrix<T> & operator-=(T const other) {
		std::for_each(body_.begin(), body_.end(),
			[&](T & value) {value -= other; });

		return *this;
	}

	//! Subtraction right matrix from left matrix
	Matrix<T> operator-(Matrix<T> const & other) const {
		// Check that rows or colls number of right and left matrix are equal
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		Matrix<T> result(*this);
	#pragma omp parallel for
		for (int i = 0; i < result.body_.size(); ++i)
			result.body_.at(i) -= other.body_.at(i);

		return result;
	}

	//! Subtraction right value from each element of left matrix
	Matrix<T> operator-(T const other) const {
		Matrix<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value -= other; });

		return result;
	}

#pragma endregion

#pragma region Operator *, *= and its overloading

	//! Multiplication each element of right matrix on value other
	Matrix<T> & operator*=(T other) {
	#pragma omp parallel for
		for (int i = 0; i < body_.size(); ++i)
			body_.at(i) *= other;

		return *this;
	}

	//! Multiplication each element of right matrix on value other
	Matrix<T> operator*(T other) const {
		Matrix<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value *= other; });

		return result;
	}

	template<typename T1>
	friend Matrix<T1> operator*(T1 const left, Matrix<T1> const & right);

	//! The scalar product of two matrices.
	//! This means that each element of left matrix multiply on corresponding element of right matrix
	Matrix<T> scalar_multiplication(Matrix<T> const & other);

#pragma endregion

#pragma region Operator /, /= and its overloading

	//! Division each element of left matrix on value other
	Matrix<T> & operator/=(T other) {
		// Check that other value is not equal by zero
		assert(other != 0);

	#pragma omp parallel for
		for (int i = 0; i < body_.size(); ++i)
			body_.at(i) /= other;

		return *this;
	}

	//! Division each element of left matrix on value other
	Matrix<T> operator/(T other) const {
		// Check that other value is not equal by zero
		assert(other != 0);

		Matrix<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value /= other; });

		return result;
	}

	//! Termwise division of the matrix argument
	//! This means each element of right matrix divide on corresponding elemrnt of left matrix
	Matrix<T> times_divide(Matrix<T> const & other);

#pragma endregion

#pragma region Properties (Get/Set methods)

	//! Get element at y row and x columns
	T & operator()(unsigned y, unsigned x) 
	{
		return body_.at(y * colls_ + x);
	}

	//! Get element at y row and x columns
	T const operator()(unsigned y, unsigned x) const 
	{
		return body_.at(y * colls_ + x);
	}

	//! Returns values pair, in witch first = rows number, second = columns number
	std::pair<unsigned int, unsigned int> size() const;

	//! Returns the sum of elements of the matrix
	long double getSum() const;

	//! Returns std::vector<T>, with values from y ID row in matrix
	std::vector<T> getRow(unsigned const y) const;

	//! Set elements of y ID row in matrix equal to values from std::vector<T> row
	void setRow(unsigned const y, std::vector<T> const & row);

	//! Returns std::vector<T>, with values in ID range [1 : rows_ - 2] column with x ID
	std::vector<T> getColumn(unsigned const x) const;

	//! Set elements x ID column of matrix in ID range [1 : colls - 1] equal to values of std::vector<T> coll
	void setColumn(unsigned const x, std::vector<T> const & coll);

#pragma endregion

#pragma region Methods

	//! Fill all matrix with value
	//! This means, that after this operation all elements of matrix are equal to value
	void fillWith(T const value);

	//! Fill elements withought boundary with value
	//! This means, that after this operation all elements besides boundaries of matrix are equal to value
	void fillWithoughtBoundary(T const value);

	//! Fill column with ID coll_id of matrix with value
	void fillColumnWith(int const coll_id, T const value);

	//! Fill row with ID row_id of matrix with value
	void fillRowWith(int const row_id, T const value);

	//! Resize matrix and fill with zeros
	void resize(unsigned rows, unsigned colls);

#pragma region Ostream Methods

	//! Print matrix with value_name name on time time step in file
	void writeToFile(std::string value_name, int const time);

	//! Print coll_id column of matrix with value_name name on time time step in file
	void writeColumnToFile(std::string value_name, int const coll_id, int const time);

	//! Print row_id row of matrix with value_name name on time time step in file
	void writeRowToFile(std::string value_name, int const row_id, int const time);

	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, Matrix<T1> const & matrix);


#pragma endregion


#pragma endregion

protected:

#pragma region Fields

	// Number of rows in matrix
	unsigned rows_;
	// Number of columns in matrix
	unsigned colls_;

	// Main body of matrix, wich contain all matrix elements
	std::vector<T> body_;

#pragma endregion

};

#pragma region additional_functions
	
#pragma endregion


#include"my_matrix_impl.h"
