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

	DistributionFunction();
	DistributionFunction(unsigned rows, unsigned colls);
	~DistributionFunction();
	
	DistributionFunction(DistributionFunction<T> const & other);


	// >> New BCs attempt 

	std::vector<T> GetTopValues(const int q) const
	{
		switch (q)
		{
		case 2:
			return dfunc_body_.at(q).GetRowMid(rows_ - 1);
			break;
		case 5:
			return dfunc_body_.at(q).GetRowRight(rows_ - 1);
			break;
		case 6:
			return dfunc_body_.at(q).GetRowLeft(rows_ - 1);
			break;
		default:
			break;
		}
	}

	// 1 - index for von-neumann bcs
	std::vector<T> GetTopValues1(const int q) const
	{
		switch (q)
		{
		case 2:
			return dfunc_body_.at(q).GetRowMid(rows_ - 2);
			break;
		case 5:
			return dfunc_body_.at(q).GetRowRight(rows_ - 2);
			break;
		case 6:
			return dfunc_body_.at(q).GetRowLeft(rows_ - 2);
			break;
			// For von neuman mid values
		case 0:
			return dfunc_body_.at(q).GetRowMid(rows_ - 2);
			break;
		case 1:
			return dfunc_body_.at(q).GetRowMid(rows_ - 2);
			break;
		case 3:
			return dfunc_body_.at(q).GetRowMid(rows_ - 2);
			break;
		default:
			break;
		}
	}

	std::vector<T> GetTopValues2(const int q) const
	{
		return dfunc_body_.at(q).GetRowMid(rows_ - 2);
	}



	std::vector<T> GetBottomValues(const int q) const
	{
		switch (q)
		{
		case 4:
			return dfunc_body_.at(q).GetRowMid(0);
			break;
		case 8:
			return dfunc_body_.at(q).GetRowRight(0);
			break;
		case 7:
			return dfunc_body_.at(q).GetRowLeft(0);
			break;
		default:
			break;
		}
	}

	std::vector<T> GetBottomValues1(const int q) const
	{
		switch (q)
		{
		case 4:
			return dfunc_body_.at(q).GetRowMid(1);
			break;
		case 8:
			return dfunc_body_.at(q).GetRowRight(1);
			break;
		case 7:
			return dfunc_body_.at(q).GetRowLeft(1);
			break;
		case 0:
			return dfunc_body_.at(q).GetRowMid(1);
			break;
		case 1:
			return dfunc_body_.at(q).GetRowMid(1);
			break;
		case 3:
			return dfunc_body_.at(q).GetRowMid(1);
			break;
		default:
			break;
		}
	}

	
	std::vector<T> GetLeftValues(const int q) const
	{
		//return dfunc_body_.at(q).GetColumn(0);

		switch (q)
		{
		case 3:
			return dfunc_body_.at(q).GetColumnMid(0);
			break;
		case 6:
			return dfunc_body_.at(q).GetColumnTop(0);
			break;
		case 7:
			return dfunc_body_.at(q).GetColumnBot(0);
			break;
		default:
			break;
		}
	}

	std::vector<T> GetLeftValues1(const int q) const
	{
		//return dfunc_body_.at(q).GetColumn(colls_ - 1);
		return dfunc_body_.at(q).GetColumnMid(1);
	}
	

	std::vector<T> GetRightValues(const int q) const
	{
		//return dfunc_body_.at(q).GetColumn(colls_ - 1);
		switch (q)
		{
		case 1:
			return dfunc_body_.at(q).GetColumnMid(colls_ - 1);
			break;
		case 5:
			return dfunc_body_.at(q).GetColumnTop(colls_ - 1);
			break;
		case 8:
			return dfunc_body_.at(q).GetColumnBot(colls_ - 1);
			break;
		default:
			break;
		}
	}

	std::vector<T> GetRightValues1(const int q) const
	{
		//return dfunc_body_.at(q).GetColumn(colls_ - 1);
		return dfunc_body_.at(q).GetColumnMid(colls_ - 2);
	}


	void SetTopValues(int const q, std::vector<T> const & row)
	{
		dfunc_body_.at(q).SetRowNew(rows_ - 2, row);
	}
	void SetBottomValues(int const q, std::vector<T> const & row)
	{
		dfunc_body_.at(q).SetRowNew(1, row);
	}
	void SetLeftValues(int const q, std::vector<T> const & coll)
	{
		dfunc_body_.at(q).SetColumnNew(1, coll, q);
	}
	void SetRightValues(int const q, std::vector<T> const & coll)
	{
		dfunc_body_.at(q).SetColumnNew(colls_ - 2, coll, q);
	}
	// >>


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