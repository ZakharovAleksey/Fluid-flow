#pragma once

#include<array>

#include"../../math/2d/my_matrix_2d.h"
#include"../../phys_values/2d/macroscopic_param_2d.h"

#include"../../modeling_area/medium.h"

// Possible number of directions in which pseudo-particles could move.
unsigned const kQ{ 9 };

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
	DistributionFunction(unsigned rows, unsigned colls);
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

	// Gets 'q'-component of probability distribution function
	Matrix2D<T> & operator[](unsigned q);
	
	// Get pair in witch: first = rows_, second = colls_
	std::pair<unsigned int, unsigned int> size() const;

	
	// —ƒ≈Ћј“№ „≈–≈« OVERRIDE „≈–≈« »Ќ“≈–‘≈…—ј

	std::vector<T> getTopBoundaryValues(int const q) const;
	std::vector<T> getBottomBoundaryValue(int const q) const;
	std::vector<T> getLeftBoundaryValue(int const q) const;
	std::vector<T> getRightBoundaryValue(int const q) const;

	
	void setTopBoundaryValue(int const q, std::vector<T> const & row);
	void setBottomBoundaryValue(int const q, std::vector<T> const & row);
	void setLeftBoundaryValue(int const q, std::vector<T> const & coll);
	void setRightBoundaryValue(int const q, std::vector<T> const & coll);

#pragma endregion

#pragma region Methods

	//! Fill each of kQ component of probability distribution function except boundaries with value
	void fillWithoutBoundaries(T const value);

	//! Fill each of kQ component of probability distribution function boundaries with value
	void fillBoundaries(T const value);

	//! Resize each of kQ component of probability distribution function 
	void resize(unsigned rows, unsigned colls);

	//! —читает плотность в кажой из €чеек области
	MacroscopicParam<T> calculateDensity() const;

	//! —читает скорость в кажой из €чеек области
	MacroscopicParam<T> calculateVelocity(const double mas[kQ], MacroscopicParam<T> const & density) const;

	//! Calculate velocity according to IB-LBM Implementation
	MacroscopicParam<T> calculateVelocity(const double mas[kQ], MacroscopicParam<T> const & density, Matrix2D<double> const & f) const;

	
	T Get(int q, int y, int x)
	{
		return dfunc_body_.at(q)(y, x);
	}

	void Set(int q, int y, int x, double value)
	{
		dfunc_body_.at(q)(y, x) = value;
	}

	void Swap(int q1, int y1, int x1, int q2, int y2, int x2)
	{
		std::swap(dfunc_body_.at(q1)(y1, x1), dfunc_body_.at(q2)(y2, x2));
	}

	const int GetRowsNumber() const;
	//! Returns number of columns
	const int GetColumnsNumber() const;


#pragma endregion

	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, DistributionFunction<T1> const & distr_func);

private:

	// Height of modeling area across Y axis direction.
	unsigned rows_;
	// Lenght of modeling area across X axis direction.
	unsigned colls_;

	// Array which store all components of probability distribution function.
	std::array<Matrix2D<T>, kQ> dfunc_body_;

};


#include"distribution_func_2d_impl.h"