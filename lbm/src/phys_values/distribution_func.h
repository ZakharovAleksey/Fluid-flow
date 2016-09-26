#pragma once

#include<array>

#include"..\math\my_matrix.h"
#include"..\phys_values\macroscopic_param.h"

//! Possible Count of directions in which particles would move
unsigned const kQ{ 9 };

/*!
	Probability distribution function implementation class.
	
	- Implemented as an array of size kQ each element of witch store probability distribution function field for apprpriate
	component of direction. 

	Example:
		DistributionFunction f[1] - stores probability distribution function field in right direction.
*/
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

	//! Get q-component of probability distribution function
	Matrix<T> & operator[](unsigned q);
	
	//! Get pair in witch: first = rows_, second = colls_
	std::pair<unsigned int, unsigned int> size() const;

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
	void resize(unsigned rows, unsigned colls);

	//! —читает плотность в кажой из €чеек области
	MacroscopicParam<T> calculateDensity() const;

	//! —читает скорость в кажой из €чеек области
	MacroscopicParam<T> calculateVelocity(const double mas[kQ], MacroscopicParam<T> const & density) const;

#pragma endregion

	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, DistributionFunction<T1> const & distr_func);

private:

#pragma region Fields
	//! Rows count for distribution function
	unsigned rows_;
	//! Columns count for distribution function
	unsigned colls_;
	
	//! Array witch store all kQ components of probability distribution function
	std::array<Matrix<T>, kQ> dfunc_body_;

#pragma endregion

};


#include"distribution_func_impl.h"