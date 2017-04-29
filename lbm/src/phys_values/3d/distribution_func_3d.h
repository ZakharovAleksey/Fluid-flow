#pragma once

#ifndef DISTRIBUTION_FUNC_3D_H
#define DISTRIBUTION_FUNC_3D_H

#include"../../math/3d/my_matrix_3d.h"
#include"../distribution_function_interface.h"

#include<array>

// Number of velocity directions in D3Q19 model
const int kQ3d = 19;


// Class: Probability distribution function in 3D case
// Consists from array of 3D matrix, which numer is equal to number of velocity directions
// in choosen model (kQ3d). 
// For example:
//	- DisributionFunc[0] - contains probability distribution function which did not move.
//  - DisributionFunc[1] - contains probability distribution function which move to the right. and so on...
//
// Changes in comparison with 2d implementation : add similar methods for new NEAR and FAR boundaries.
template<typename T>
class DistributionFunction3D : public iDistributionFunction<T>
{
public:
	DistributionFunction3D(int depth, int rows, int colls);
	~DistributionFunction3D();

#pragma region Operators overloading

	DistributionFunction3D<T> & operator=(const DistributionFunction3D<T> & other)
	{
		if (this != &other)
		{
			DistributionFunction3D<T> res(other);
			Swap(res);
		}

		return *this;
	}

#pragma region Unary operators

	template<typename T1>
	friend const DistributionFunction3D<T1> operator-(const DistributionFunction3D<T> & right);

#pragma endregion

#pragma region Binary operators

	template<typename T1>
	friend DistributionFunction3D<T1> & operator+=(DistributionFunction3D<T1> & left, const DistributionFunction3D<T1> & right);

	template<typename T1>
	friend DistributionFunction3D<T1> & operator+=(DistributionFunction3D<T1> & left, const T1 & right);

	template<typename T1>
	friend DistributionFunction3D<T1> & operator-=(DistributionFunction3D<T1> & left, const DistributionFunction3D<T1> & right);

	template<typename T1>
	friend DistributionFunction3D<T1> & operator-=(DistributionFunction3D<T1> & left, const T1 & right);
	
	template<typename T1>
	friend DistributionFunction3D<T1> & operator*=(DistributionFunction3D<T1> & left, const T1 & right);

	template<typename T1>
	friend DistributionFunction3D<T1> & operator/=(DistributionFunction3D<T1> & left, const T1 & right);

#pragma endregion

#pragma region Streaming operators

	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, const DistributionFunction3D<T1> & dist_func);

#pragma endregion

	// Gets 'q'-s component of distribution function (As reference)
	Matrix3D<T> & operator[](const int q);
	// Gets 'q'-s component of distribution function (As constant value)
	const Matrix3D<T> & operator[](const int q) const;

#pragma endregion


	// Swaps current matrix with 'other' matrix
	void Swap(DistributionFunction3D<T> & other);

	// Override methods of iDistributionFunction 
	std::vector<T> GetTopBoundaryValues(int const q) const  override;
	std::vector<T> GetBottomBoundaryValue(int const q) const override;
	std::vector<T> GetLeftBoundaryValue(int const q) const override;
	std::vector<T> GetRightBoundaryValue(int const q) const override;

	//! Gets probability distribution function values on NEAR boundary
	std::vector<T> GetNearBoundaryValue(int const q) const;
	//! Gets probability distribution function values on FAR boundary
	std::vector<T> GetFarBoundaryValue(int const q) const;

	// Override methods of iDistributionFunction 
	void SetTopBoundaryValue(int const q, std::vector<T> const & layer) override;
	void SetBottomBoundaryValue(int const q, std::vector<T> const & layer) override;
	void SetLeftBoundaryValue(int const q, std::vector<T> const & layer) override;
	void SetRightBoundaryValue(int const q, std::vector<T> const & layer) override;

	//! Set probability distribution function values on NEAR boundary equal to 'layer'
	void SetNearBoundaryValue(int const q, std::vector<T> const & layer);
	//! Set probability distribution function values on FAR boundary equal to 'layer'
	void SetFarBoundaryValue(int const q, std::vector<T> const & layer);

	//! Clear all boundaries of distribution function
	void ClearBoundaries();


private:
	//! Depth of the matrix (Z-axis size  value)
	int depth_;
	//! Number of matrix rows (Y-axis size  value)
	int rows_;
	//! Number of matrix columns (X-axis size  value)
	int colls_;

	//! Main body of distribution function, containing all it's components
	std::array<Matrix3D<T>, kQ3d> body_;
};


#include"distribution_func_3d_impl.h"

#endif // !DISTRIBUTION_FUNC_3D_H

