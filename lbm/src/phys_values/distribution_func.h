#pragma once

#include<array>

#include"..\math\my_matrix_2d.h"
#include"..\phys_values\macroscopic_param.h"

/// <summary>
/// Possible number of directions in which pseudo-particles could move.
/// </summary>
int const kQ{ 9 };

/// <summary>
/// Template class for probability distribution function field for all modeling area representation.
/// <remarks> 
/// Let us assume, we deal with velocity model D2Q9 for Lattice Boltzmann method.
/// 
/// This means we have 2D model with 9 directions in which pseudo-particles could move.
/// So, for this velocity model current class field represent from itself an array of 9 
/// elements. 
/// Each of this 9 elements represent from itself appropriate component of probability 
/// distribution function field (size of each is equal to modeling area size).
/// This means:
///  - body[0] - store 0-distribution function component field.
///  - dfunc_body[1] - store 1-distribution function component field.
///  ... and so on.
/// </remarks>
/// <example>
/// blah blah
/// </example>
/// </summary>
template<typename T>
class DistributionFunction
{
public:

#pragma region Constructor

	DistributionFunction();
	DistributionFunction(int rows, int colls);
	~DistributionFunction();
	
	DistributionFunction(DistributionFunction<T> const & other);

#pragma endregion

#pragma region Overload operators

	void swap(DistributionFunction & dist_func);

	DistributionFunction<T> & operator=(DistributionFunction<T> const & other) {
		assert(rows_ == other.rows_ && colls_ == other.colls_);
		if (this != &other) {
			DistributionFunction<T> temp(other);
			temp.swap(*this);
		}
		return *this;
	}

#pragma region Operator +, += overload

	DistributionFunction<T> & operator+=(DistributionFunction<T> const & other) 
	{
		for (int q = 0; q < kQ; ++q)
			dfunc_body_.at(q) += other.dfunc_body_.at(q);

		return *this;
	}

	DistributionFunction<T> operator+(DistributionFunction<T> const & other) 
	{
		DistributionFunction<T> result(*this);
		for (int q = 0; q < kQ; ++q)
			result.dfunc_body_.at(q) += other.dfunc_body_.at(q);

		return result;
	}

#pragma endregion

#pragma region  Operator -, -= overload

	DistributionFunction<T> & operator-=(DistributionFunction<T> const & other) 
	{
		for (int q = 0; q < kQ; ++q)
			dfunc_body_.at(q) -= other.dfunc_body_.at(q);

		return *this;
	}

	DistributionFunction<T> operator-(DistributionFunction<T> const & other) 
	{
		DistributionFunction<T> result(*this);

		for (int q = 0; q < kQ; ++q)
			result.dfunc_body_.at(q) -= other.dfunc_body_.at(q);

		return result;
	}

#pragma endregion

#pragma region Operator / , /= overload
	
	DistributionFunction<T> & operator/=(T const & other) 
	{
		for (int q = 0; q < kQ; ++q)
			dfunc_body_.at(q) /= other;

		return *this;
	}

	DistributionFunction<T> operator/(T const & other) 
	{
		DistributionFunction<T> result(*this);

		for (int q = 0; q < kQ; ++q)
			result.dfunc_body_.at(q) /= other;

		return result;
	}

#pragma endregion

#pragma endregion

#pragma region Proprerties (Get/Set)

	//! Get q-component of probability distribution function
	Matrix2D<T> & operator[](int q);
	
	//! Get pair in witch: first = rows_, second = colls_
	std::pair<int, int> size() const;

	//! Get probability distribution function values on TOP boundary
	std::vector<T> getTopBoundaryValues(int const q) const;
	//! Get probability distribution function values on BOTTOM boundary
	std::vector<T> getBottomBoundaryValue(int const q) const;
	//! Get probability distribution function values on LEFT boundary
	std::vector<T> getLeftBoundaryValue(int const q) const;
	//! Get probability distribution function values on RIGHT boundary
	std::vector<T> getRightBoundaryValue(int const q) const;

	//! Set probability distribution function values on TOP boundary equal to parameter array
	void setTopBoundaryValue(int const q, std::vector<T> const & row);
	//! Set probability distribution function values on BOTTOM boundary equal to parameter array
	void setBottomBoundaryValue(int const q, std::vector<T> const & row);
	//! Set probability distribution function values on LEFT boundary equal to parameter array
	void setLeftBoundaryValue(int const q, std::vector<T> const & coll);
	//! Set probability distribution function values on RIGHT boundary equal to parameter array
	void setRightBoundaryValue(int const q, std::vector<T> const & coll);

#pragma endregion

#pragma region Methods

	//! Fill each of kQ component of probability distribution function except boundaries with value
	void fillWithoutBoundaries(T const value);

	//! Fill each of kQ component of probability distribution function boundaries with value
	void fillBoundaries(T const value);

	//! Resize each of kQ component of probability distribution function 
	void resize(int rows, int colls);

	//! —читает плотность в кажой из €чеек области
	MacroscopicParam<T> calculateDensity() const;

	//! —читает скорость в кажой из €чеек области
	MacroscopicParam<T> calculateVelocity(const double mas[kQ], MacroscopicParam<T> const & density) const;

#pragma endregion

	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, DistributionFunction<T1> const & distr_func);

private:

#pragma region Fields

	/// <summary>
	/// Height of modeling area across Y axis direction.
	/// </summary>
	int rows_;
	
	/// <summary>
	/// Lenght of modeling area across X axis direction.
	/// </summary>
	int colls_;

	/// <summary>
	/// Array which store all components of probability distribution function.
	/// </summary>
	std::array<Matrix2D<T>, kQ> dfunc_body_;

#pragma endregion

};


#include"distribution_func_impl.h"