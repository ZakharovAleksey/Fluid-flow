#pragma once

#ifndef MY_MATRIX_3D_H
#define MY_MATRIX_3D_H

#include"../my_matrix_interface.h"
#include"../2d/my_matrix_2d.h"

template<typename T>
class Matrix3D : public iMatrix<T>
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


	// Get Top or Bottom layer for 'z' is equal to 0 or depth-1
	std::vector<T> GetTBLayer(const int z) const;
	// Get Left or Right layer for 'x' is equal to 0 or colls-1 (without upper and down elements)
	std::vector<T> GetLRLayer(const int x) const;
	// Get Left or Right layer for 'y' is equal to 0 or rows-1 (without upper and down elements)
	std::vector<T> GetNFLayer(const int y) const;

	// Set Top or Bottom layer for 'z' is equal to 0 or depth-1 to 'layer'
	void SetTBLayer(unsigned const z, std::vector<T> const & layer);
	// Set Left or Right layer for 'x' is equal to 0 or colls-1 to 'layer' (without upper and down elements)
	void SetLRLayer(unsigned const x, std::vector<T> const & layer);
	// Set Left or Right layer for 'y' is equal to 0 or rows-1 to 'layer' (without upper and down elements)
	void SetNFLayer(unsigned const y, std::vector<T> const & layer);

	// Override of iMatrix method
	long double GetSum() const override;


#pragma endregion

	//! Scalar multiplication on 'other' matrix: each element of matrix multiplies on appropriate element
	// of 'other' matrix : this[1][2][3] * other[1][2][3]
	Matrix3D<T> const ScalarMultiplication(Matrix3D<T> const & other);
	//! Times divide on 'other' matrix: each element of matrix divides on appropriate element
	// of 'other' matrix : this[1][2][3] / other[1][2][3]
	Matrix3D<T> const TimesDivide(Matrix3D<T> const & other);

	
	//! Override of iMatrix method [As far as I rememder did not used in code because of std::unique_ptr<>]
	void Resize(int new_rows_numb, int new_colls_numb, int new_depth_numb = 0) override;

	//! Filles ONLY side walls of Matrix with 'value' (not TOP and BOTTOM)
	void FillBoundarySideWalls(const T value);
	//! Filles ONLY ONE layer of Matrix with 'value' (All Oxy plane with fixed 'z')
	void FillLayer(const int z, const T value);
	//! Filles ONLY TOP and BOTTOM layers with 'value'
	void FillTopBottomWalls(const T value);
	//! Filles all elements of matrix besides it's boundaries with 'value'
	void FillWithoutBoundary(const T value);
	//! Filles all elements of matrix with 'value'
	void FillWith(const T value);


private:

	//! Swaps current 3d-matrix with 'other' 3d-matrix (No any check on dimension sizes)
	void Swap(Matrix3D<T> & other);
	//! Returns total number of elements in 3D matrix
	int const GetTotalSize() const { return depth_ * rows_ * colls_; }

private:
	//! Depth of the matrix (Z-axis size  value)
	int depth_;
	//! Number of matrix rows (Y-axis size  value)
	int rows_;
	//! Number of matrix columns (X-axis size  value)
	int colls_;

	//! Body of matrix, stores all matrix elements
	std::vector<T> body_;

};


#include"my_matrix_3d_impl.h"

#endif // !MY_MATRIX_3D_H