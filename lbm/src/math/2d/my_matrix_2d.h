#pragma once

#include"../my_matrix_interface.h"

#include<iostream>
#include<ostream>
#include <cassert> // setw()
#include<vector>
#include <functional>
#include<string>

#include<ppl.h>
#include<omp.h>


#define TEST_INCLUDE 

#ifndef TEST_INCLUDE 
 	#include<cstdlib>
 	#include<ctime>
#endif // !TEST_INCLUDE


/// <summary>
/// The template class of 2D matrix.
/// </summary>
/// <remarks>
/// Basic class for LBM project.
/// All physical values (density, velocity, distribution function) will be implemented as derived clases from Matrx class.
/// Overloading implemetation only for operations witch we need in LBM implementation.
/// </remarks>
template<typename T>
class Matrix2D : public iMatrix<T>
{
public:

#pragma region Constructor

	/// <summary>
	/// Allocates an uninitialized matrix that holds zero elements.
	/// </summary>
	Matrix2D();

	/// <summary>
	/// Allocates an matrix with rows number is equal to "rows" and column number is equal to "cells", which is filled with zeros.
	/// </summary>
	/// <param name="rows"> Number of rows in matrix. </param>
	/// <param name="colls"> Number of columns in matrix.</param>
	Matrix2D(unsigned rows, unsigned colls);

	virtual ~Matrix2D();

	/// <summary>
	/// Allocates an matrix with rows number is equal to the "other" matrix rows and column number is equal to the "other" matrix cells, which is filled with 
	/// the "other" matrix values.
	/// </summary>
	Matrix2D(Matrix2D<T> const & other);	// ��������� ��� �� ������ ��� ������ ������� ������

	/// <summary>
	/// Swap current matrix (Swap includes rows number, columns number and all elements) with the "other" matrix.
	/// </summary>
	/// <param name="other"> The matrix which is swapped with the current </param>
	void Swap(Matrix2D<T> & other);

#pragma endregion

	/// <summary>
	/// The assignment operator.
	/// </summary>
	/// <param name="other"> Matirix from the rigth side of the assigment operator. </param>
	/// <returns> New assigment matrix. </returns>
	Matrix2D<T> & operator=(Matrix2D<T> const & other) 
	{
		// Check that rows or colls number of right and left matrix are equal
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		if (this != & other) 
		{
			Matrix2D<T> temp(other);
			temp.Swap(*this);
		}

		return *this;
	}

#pragma region Operator +, += and its overloading

	/// <summary>
	/// A computed assignment operator. Adds the matrix "other" expression to the matrix.
	/// </summary>
	/// <param name="other"> Matrix that we add to the current matrix. </param>
	/// <returns>  The result of adding two matrices. </returns>
	Matrix2D<T> & operator+=(Matrix2D<T> const & other) 
	{
		assert(rows_ == other.rows_ && colls_ == other.colls_);
	
		#pragma omp parallel for
			for (int i = 0; i < body_.size(); ++i)
				body_.at(i) += other.body_.at(i);

		return *this;
	}

	/// <summary>
	/// A computed assignment operator. Adds the to each matrix element "other" value.
	/// </summary>
	/// <param name="other"> Value, we add to the each element of the current matrix. </param>
	/// <returns>  The result of adding the matrix and the value.  </returns>
	Matrix2D<T> & operator+=(T const other) 
	{
		std::for_each(body_.begin(), body_.end(), [&](T & value) {value += other; });
		return *this;
	}

	/// <summary>
	/// A computed assignment operator. Adds the current matrix with the "other" matrix.
	/// </summary>
	/// <param name="other"> Matrix that we add to the current matrix. </param>
	/// <returns>  The result of adding two matrices. </returns>
	Matrix2D<T> operator+(Matrix2D<T> const & other) const 
	{
		// Check that rows or columns number of right and left matrix are equal
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		Matrix2D<T> result(*this);
	#pragma omp parallel for
		for (int i = 0; i < result.body_.size(); ++i)
			result.body_.at(i) += other.body_.at(i);

		return result;
	}

	/// <summary>
	/// A computed assignment operator. Adds the to each element of the current matrix element "other" value.
	/// </summary>
	/// <param name="other"> Value, we add to the each element of the current matrix. </param>
	/// <returns> The result of adding the matrix and the value. </returns>
	Matrix2D<T> operator+(T const other) const 
	{
		Matrix2D<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(), [&](T & value) {value += other; });

		return result;
	}

	
	template<typename T1>
	/// <summary>
	/// A computed assignment operator. Adds the to each element of the current matrix element "other" value.
	/// </summary>
	/// <param name="left"> Value, we add to the each element of the current matrix. </param>
	/// <param name="right"> the current matrix.</param>
	/// <returns> The result of adding the value and the matrix. </returns>
	friend Matrix2D<T1> operator+(T1 const left, Matrix2D<T1> const & right);

#pragma endregion

#pragma region Operator -, -= and its overloading

	/// <summary>
	/// A computed assignment operator. Substract the matrix "other" expression to the matrix.
	/// </summary>
	/// <param name="other"> Matrix that we substract from the current matrix. </param>
	/// <returns>  The result of subctract two matrices. </returns>
	Matrix2D<T> & operator-=(Matrix2D<T> const & other) {
		// Check that rows or columns number of right and left matrix are equal
		assert(rows_ == other.rows_ && colls_ == other.colls_);
	#pragma omp parallel for
			for (int i = 0; i < body_.size(); ++i)
				body_.at(i) -= other.body_.at(i);

		return *this;
	}

	/// <summary>
	/// A computed assignment operator. Subctract the to each matrix element "other" value.
	/// </summary>
	/// <param name="other"> Value, we substract from the each element of the current matrix. </param>
	/// <returns>  The result of subctracting the matrix and the value.  </returns>
	Matrix2D<T> & operator-=(T const other) {
		std::for_each(body_.begin(), body_.end(),
			[&](T & value) {value -= other; });

		return *this;
	}

	/// <summary>
	/// A computed assignment operator. Substract the "other" matrix from the current matrix.
	/// </summary>
	/// <param name="other"> Matrix that we substract from the current matrix. </param>
	/// <returns>  The result of substracting two matrices. </returns>
	Matrix2D<T> operator-(Matrix2D<T> const & other) const {
		// Check that rows or colls number of right and left matrix are equal
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		Matrix2D<T> result(*this);
	#pragma omp parallel for
		for (int i = 0; i < result.body_.size(); ++i)
			result.body_.at(i) -= other.body_.at(i);

		return result;
	}

	/// <summary>
	/// A computed assignment operator. Substract from the each element of the current matrix element "other" value.
	/// </summary>
	/// <param name="other"> Value, we substract from the each element of the current matrix. </param>
	/// <returns> The result of adding the matrix and the value. </returns>
	Matrix2D<T> operator-(T const other) const {
		Matrix2D<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value -= other; });

		return result;
	}

#pragma endregion

#pragma region Operator *, *= and its overloading

	/// <summary>
	/// A computed assignment operator. Multiplies each element of the matrix with "other" value.
	/// </summary>
	/// <param name="other"> Value that we multiply on the current matrix. </param>
	/// <returns>  The result of multiplying matrix with value. </returns>
	Matrix2D<T> & operator*=(T other) {
	#pragma omp parallel for
		for (int i = 0; i < body_.size(); ++i)
			body_.at(i) *= other;

		return *this;
	}

	/// <summary>
	///  A computed assignment operator. Multiplies each element of the matrix with "other" value.
	/// </summary>
	/// <param name="other"> Value that we multiply to the current matrix. </param>
	/// <returns>  The result of multiplying the matrix with value. </returns>
	Matrix2D<T> operator*(T other) const {
		Matrix2D<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value *= other; });

		return result;
	}

	template<typename T1>
	/// <summary>
	/// A computed assignment operator. Multiplies each element of the matrix "other" expression with "other" value.
	/// </summary>
	/// <param name="left"> Value that we multiply to the current matrix. </param>
	/// <param name="right"> The current matrix. </param>
	/// <returns> The result of multiplying the matrix with value. </returns>
	friend Matrix2D<T1> operator*(T1 const left, Matrix2D<T1> const & right);

	/// <summary>
	/// The scalar product of current matrix and "other" matrix.
	/// </summary>
	/// <param name="other"> Right part of scalar product matrix. </param>
	/// <returns> The scalar product of two matrix. </returns>
	Matrix2D<T> ScalarMultiplication(Matrix2D<T> const & other);

#pragma endregion

#pragma region Operator /, /= and its overloading

	/// <summary>
	/// A computed assignment operator. Divide each element of the matrix on "other" value.
	/// </summary>
	/// <param name="other"> Matrix that we divide on the current matrix. </param>
	/// <returns>  The result of divideing two matrices. </returns>
	Matrix2D<T> & operator/=(T other) {
		// Check that other value is not equal by zero
		assert(other != 0);

	#pragma omp parallel for
		for (int i = 0; i < body_.size(); ++i)
			body_.at(i) /= other;

		return *this;
	}

	/// <summary>
	/// A computed assignment operator. Divide each element of the matrix on "other" value.
	/// </summary>
	/// <param name="other"> Value by which we divide each element of the matrix. </param>
	/// <returns>  The result of divideing two matrices. </returns>
	Matrix2D<T> operator/(T other) const {
		// Check that other value is not equal by zero
		assert(other != 0);

		Matrix2D<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value /= other; });

		return result;
	}

	/// <summary>
	/// Termwise division of the matrix on "other" argument.
	/// </summary>
	/// <param name="other"> Matrix by which we divide aproppriate element of the current matrix. </param>
	/// <returns> Termwise division of the two matrix. </returns>
	Matrix2D<T> TimesDivide(Matrix2D<T> const & other);

#pragma endregion

#pragma region Properties (Get/Set methods)

	/// <summary>
	/// Get element at choosen row and column respectively
	/// </summary>
	/// <param name="y"> Row index </param>
	/// <param name="x"> Column index </param>
	/// <returns> Value in matrix at [y][x] position</returns>
	T & operator()(unsigned y, unsigned x) 
	{
		return body_.at(y * colls_ + x);
	}

	/// <summary>
	/// Get element at choosen row and column respectively as a constant
	/// </summary>
	/// <param name="y"> Row index </param>
	/// <param name="x"> Column index </param>
	/// <returns></returns>
	T const operator()(unsigned y, unsigned x) const 
	{
		return body_.at(y * colls_ + x);
	}

	

	/// <summary>
	/// Returns values pair, in witch first isequal to rows number, second is equal to columns number
	/// </summary>
	/// <returns> Pair, where first element is row number and second is column number </returns>
	std::pair<unsigned int, unsigned int> Size() const;

	/// <summary>
	/// Returns the sum of all elements of the matrix
	/// </summary>
	/// <returns></returns>
	long double GetSum() const;

	/// <summary>
	/// Returns vector, with values from row with choosen index from matrix
	/// </summary>
	/// <param name="y"> Index of row </param>
	/// <returns> vector with values from choosen row </returns>
	std::vector<T> GetRow(unsigned const y) const;

	/// <summary>
	/// Set elements of "y" row in the current matrix equal to values from std::vector "row"
	/// </summary>
	/// <param name="y"> Index of matrix row for insertion operation. </param>
	/// <param name="row"> Vector for for insertion operation. </param>
	void SetRow(unsigned const y, std::vector<T> const & row);

	/// <summary>
	/// Returns std::vector, with values in range [1 : rows_ - 2] from column with "x" index.
	/// </summary>
	/// <param name="x"> Index of column from witch method extract values. </param>
	/// <returns> Vector filled with values from apropriate column of the current matrix. </returns>
	std::vector<T> GetColumn(unsigned const x) const;

	/// <summary>
	/// Set elements of "x" column in the current matrix equal to values from std::vector "coll"
	/// </summary>
	/// <param name="x"> Index of matrix column for insertion operation. </param>
	/// <param name="coll">  Vector for for insertion operation. </param>
	void SetColumn(unsigned const x, std::vector<T> const & coll);

#pragma endregion

#pragma region Methods

	/// <summary>
	/// Fills the matrix values with "value" parameter.
	/// </summary>
	/// <param name="value"> The value that is filled elements of the the matrix. </param>
	void FillWith(T const value);

	/// <summary>
	/// Fills the matrix values withought boundary with "value" parameter. 
	/// This means, that after this operation all elements besides boundaries of the matrix are equal to "value" parameter.
	/// </summary>
	/// 
	/// <param name="value"> The value that is filled elements of the the matrix. </param>
	void FillWithoughtBoundary(T const value);

	/// <summary>
	/// Fill "coll_id" column in the matrix with "value" parameter.
	/// </summary>
	/// <param name="coll_id"> �olumn index in the matrix to be filled with "value" parameter. </param>
	/// <param name="value"> The value that is filled elements of the the matrix. </param>
	void FillColumnWith(int const coll_id, T const value);

	/// <summary>
	/// Fill "row_id" row in the matrix with "value" parameter.
	/// </summary>
	/// <param name="row_id"> Row index in the matrix to be filled with "value" parameter. </param>
	/// <param name="value"> The value that is filled elements of the the matrix. </param>
	void FillRowWith(int const row_id, T const value);


	void Resize(int new_rows_numb, int new_colls_numb, int new_depth_numb = 0) override;

#pragma region Ostream Methods

	/// <summary>
	/// Write down selected macroscopic physical value "value_type" at iteration number "time" to *.txt file.
	/// </summary>
	/// <param name="value_name"> Selected macroscopic physical value (density or velocity). </param>
	/// <param name="time"> Selected iteration number. </param>
	void WriteToFile(std::string value_name, int const time);

	void WriteFieldToTxt(std::string path, std::string phys_val, const int time)
	{
		std::string file_name = path + "\\" + phys_val + "[" + std::to_string(rows_) + "x" + std::to_string(colls_) + "]_t" + std::to_string(time) + ".txt";

		std::ofstream output_file;
		output_file.open(file_name);

		if (output_file.is_open())
		{
			unsigned cur_position{ 1 };

			for (unsigned y = 0; y < rows_; ++y)
			{
				for (unsigned x = 0; x < colls_; ++x)
				{
					output_file << this->operator()(y, x);
					if (x != colls_ - 1)
						output_file << ' ';
				}
				if (y != rows_ - 1)
					output_file << std::endl;

			}
		}
		else
		{
			std::cout << "Error! Could not open file " << file_name << " to write data.\n";
		}

		output_file.close();
	}

	/// <summary>
	/// Write down "coll_id" column of selected macroscopic physical value "value_type" at iteration number "time" to *.txt file.
	/// </summary>
	/// <param name="value_name">  Selected macroscopic physical value (density or velocity). </param>
	/// <param name="coll_id"> The selected column index. </param>
	/// <param name="time"> Selected iteration number. </param>
	void WriteColumnToFile(std::string value_name, int const coll_id, int const time);

	/// <summary>
	/// Write down "coll_id" row of selected macroscopic physical value "value_type" at iteration number "time" to *.txt file.
	/// </summary>
	/// <param name="value_name">  Selected macroscopic physical value (density or velocity). </param>
	/// <param name="row_id"> The selected column index. </param>
	/// <param name="time"> Selected iteration number. </param>
	void WriteRowToFile(std::string value_name, int const row_id, int const time);

	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, Matrix2D<T1> const & matrix);


#pragma endregion


#pragma endregion

protected:

	// Number of rows in matrix
	int rows_;
	// Number of columns in matrix
	int colls_;
	// Main body of matrix, wich contain all matrix elements
	std::vector<T> body_;

};



#include"my_matrix_2d_impl.h"
