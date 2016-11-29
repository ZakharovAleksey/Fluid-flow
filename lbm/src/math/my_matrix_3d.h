#pragma once

#include<iostream>
#include<vector>

#include"my_matrix_interface.h"
#include"my_matrix_2d.h"


template<typename T>
class Matrix3D : public Matrix2D<T>
{
public:

	Matrix3D(const int rows, const int colls, const int height);
	~Matrix3D() {}

#pragma region Overloading functions

	/// <summary>
	/// Returns the sum of all elements of the matrix
	/// </summary>
	/// <returns></returns>
	long double GetSum() const { return 1.0; }

	/// <summary>
	/// Returns vector, with values from row with choosen index from matrix
	/// </summary>
	/// <param name="y"> Index of row </param>
	/// <returns> vector with values from choosen row </returns>
	std::vector<T> GetRow(unsigned const y) const { std::vector<int> test; return test; }

	/// <summary>
	/// Set elements of "y" row in the current matrix equal to values from std::vector "row"
	/// </summary>
	/// <param name="y"> Index of matrix row for insertion operation. </param>
	/// <param name="row"> Vector for for insertion operation. </param>
	void SetRow(unsigned const y, std::vector<T> const & row) {}

	/// <summary>
	/// Returns std::vector, with values in range [1 : rows_ - 2] from column with "x" index.
	/// </summary>
	/// <param name="x"> Index of column from witch method extract values. </param>
	/// <returns> Vector filled with values from apropriate column of the current matrix. </returns>
	std::vector<T> GetColumn(unsigned const x) const { std::vector<int> test; return test; }

	/// <summary>
	/// Set elements of "x" column in the current matrix equal to values from std::vector "coll"
	/// </summary>
	/// <param name="x"> Index of matrix column for insertion operation. </param>
	/// <param name="coll">  Vector for for insertion operation. </param>
	void SetColumn(unsigned const x, std::vector<T> const & coll) {}

#pragma endregion


	//T const operator()(const int x, const int y, const int z) const {}

	//T const Get(const int x, const int y, const int z) const;

	//T & operator()(const int x, const int y, const int z) { return T(); }

	void Set(const int x, const int y, const int z, T const value);


	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, Matrix3D<T1> const & matrix);


#pragma region Properties

#pragma endregion



private:

	int rows_;
	int colls_;
	int height_;

	typedef std::vector<Matrix2D<T>> VectorOf2DMatrix;

	VectorOf2DMatrix body3d_;

};

#include"my_matrix_3d_impl.h"
