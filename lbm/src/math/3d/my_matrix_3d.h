#pragma once

#ifndef MY_MATRIX_3D_H
#define MY_MATRIX_3D_H

#include"../my_matrix_interface.h"
#include"../2d/my_matrix_2d.h"

template<typename T>
class Matrix3D : public Matrix2D<T>, public iMatrix<T>
{
public:

	Matrix3D();
	Matrix3D(int rows, int colls, int depth);

	~Matrix3D();


#pragma region Operators overloading
	
	Matrix3D<T> & operator=(const Matrix3D<T> & other)
	{

		if (this != &other)
		{
			CompareSize(*this, other);
			Matrix3D<T> temp(other);
			Swap(temp);
		}

		return *this;
	}
	
#pragma region Unary operators

	template<typename T1>
	friend const Matrix3D<T1> operator-(const Matrix3D<T1> & right);

#pragma endregion


#pragma region Binary operators

	template<typename T1>
	friend Matrix3D<T1> & operator+=(Matrix3D<T1> & left, const Matrix3D<T1> & right);

	template<typename T1>
	friend Matrix3D<T1> & operator+=(Matrix3D<T1> & left, const T1 & value);

	template<typename T1>
	friend Matrix3D<T1> & operator-=(Matrix3D<T1> & left, const Matrix3D<T1> & right);

	template<typename T1>
	friend Matrix3D<T1> & operator-=(Matrix3D<T1> & left, const T1 & right);


	template<typename T1>
	friend Matrix3D<T1> & operator*=(Matrix3D<T1> & left, const T1 & right);

	template<typename T1>
	friend Matrix3D<T1> & operator/=(Matrix3D<T1> & left, const T1 & right);

#pragma endregion


#pragma region Stream operators

	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, const Matrix3D<T1> & matrix);

#pragma endregion


#pragma endregion
	
#pragma region Getters and Setters

	// Gets value of matrix, which correspond to the input z,y and x values (As reference)
	T & operator()(int z, int y, int x)
	{
		return body_.at(z * rows_ * colls_ + y * colls_ + x);
	}
	// Gets value of matrix, which correspond to the input z,y and x values (As constant)
	T const operator()(int z, int y, int x) const
	{
		return body_.at(z * rows_ * colls_ + y * colls_ + x);
	}

	// Returns Z-dimension size of the current 3d-matrix
	int const GetDepthNumber() const { return depth_; }
	// Returns Y-dimension size of the current 3d-matrix
	int const GetRowsNumber() const { return rows_; }
	// Returns X-dimension size of the current 3d-matrix
	int const GetCollsNumber() const { return colls_; }

#pragma endregion

	Matrix3D<T> const ScalarMultiplication(Matrix3D<T> const & other);



	long double GetSum() const override;

	std::vector<T> GetRow(unsigned const y) const override;

	void SetRow(unsigned const y, std::vector<T> const & row) override;

	std::vector<T> GetColumn(unsigned const x) const override;

	// Set elements of "x" column in the current matrix equal to values from std::vector "coll"
	void SetColumn(unsigned const x, std::vector<T> const & coll) override;





private:

	// Swaps current 3d-matrix with 'other' 3d-matrix (No any check on dimension sizes)
	void Swap(Matrix3D<T> & other);
	// Returns total number of elements in 3D matrix
	int const GetTotalSize() const { return depth_ * rows_ * colls_; }

private:

	// Depth of the matrix (Z-axis size  value)
	int depth_;
};


#include"my_matrix_3d_impl.h"

#endif // !MY_MATRIX_3D_H