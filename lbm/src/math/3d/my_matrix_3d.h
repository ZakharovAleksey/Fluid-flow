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

#pragma region Override interface methods

	long double GetSum() const override;

	std::vector<T> GetRow(unsigned const y) const override;
	std::vector<T> GetColumn(unsigned const x) const override;

	void SetRow(unsigned const y, std::vector<T> const & row) override;
	void SetColumn(unsigned const x, std::vector<T> const & coll) override;

#pragma endregion

	// Get Top or Bottom layer
	std::vector<T> GetTBLayer(const int z) const;
	// Get Left or Right layer (without upper and down elements)
	std::vector<T> GetLRLayer(const int x) const;
	// Get Left or Right layer (without upper and down elements)
	std::vector<T> GetNFLayer(const int y) const;

	void SetTBLayer(unsigned const z, std::vector<T> const & layer);
	void SetLRLayer(unsigned const x, std::vector<T> const & layer);
	void SetNFLayer(unsigned const y, std::vector<T> const & layer);




#pragma endregion

	Matrix3D<T> const ScalarMultiplication(Matrix3D<T> const & other);
	Matrix3D<T> const TimesDivide(Matrix3D<T> const & other);

	
	void Resize(int new_rows_numb, int new_colls_numb, int new_depth_numb = 0) override;

	// Filles ONLY side walls of Matrix with 'value' (not TOP and BOTTOM)
	void FillBoundarySideWalls(T value)
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

	// Filles ONLY ONE layer of Matrix with 'value' (All Oxy plane with fixed 'z')
	void FillLayer(const int z, const T value)
	{
		for (int y = 0; y < rows_; ++y)
			for (int x = 0; x < colls_; ++x)
			{
				this->operator()(z, y, x) = value;
			}
	}

	// Filles ONLY TOP and BOTTOM layers with 'value'
	void FillTopBottomWalls(const T value)
	{
		FillLayer(0, value);
		FillLayer(depth_ - 1, value);
	}


	//  !!! Чисто для тестов - потом убратьЫ !!!
	void FillWithoutBoundary(T value)
	{
		for (int z = 1; z < depth_ - 1; ++z)
			for (int y = 1; y < rows_ - 1; ++y)
				for (int x = 1; x < colls_ - 1; ++x)
					body_.at(z * rows_ * colls_ + y * colls_ + x) = value;
		
	}
	// !!! 
	void FillWith(T value)
	{
		for (auto & i : body_)
			i = value;
	}


private:

	// Swaps current 3d-matrix with 'other' 3d-matrix (No any check on dimension sizes)
	void Swap(Matrix3D<T> & other);
	// Returns total number of elements in 3D matrix
	int const GetTotalSize() const { return depth_ * rows_ * colls_; }

private:
	// Depth of the matrix (Z-axis size  value)
	int depth_;
	int rows_;
	int colls_;

	std::vector<T> body_;

};


#include"my_matrix_3d_impl.h"

#endif // !MY_MATRIX_3D_H