#pragma once

#ifndef DISTRIBUTION_FUNC_3D_H
#define DISTRIBUTION_FUNC_3D_H

#include"../../math/3d/my_matrix_3d.h"
#include"../distribution_function_interface.h"

#include<array>

const int kQ3d = 2;



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


	std::vector<T> GetTopBoundaryValues(int const q) const  override;
	std::vector<T> GetBottomBoundaryValue(int const q) const override;
	std::vector<T> GetLeftBoundaryValue(int const q) const override;
	std::vector<T> GetRightBoundaryValue(int const q) const override;

	
	void SetTopBoundaryValue(int const q, std::vector<T> const & row) override;
	void SetBottomBoundaryValue(int const q, std::vector<T> const & row) override;
	void SetLeftBoundaryValue(int const q, std::vector<T> const & coll) override;
	void SetRightBoundaryValue(int const q, std::vector<T> const & coll) override;


private:
	int depth_;
	int rows_;
	int colls_;

	std::array<Matrix3D<T>, kQ3d> body_;
};


#include"distribution_func_3d_impl.h"

#endif // !DISTRIBUTION_FUNC_3D_H

