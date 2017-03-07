#pragma once

#ifndef DISTRIBUTION_FUNC_INTERFACE_H
#define DISTRIBUTION_FUNC_INTERFACE_H

template<typename T>
class iDistributionFunction
{
public:
	virtual ~iDistributionFunction() {}


	virtual std::vector<T> getTopBoundaryValues(int const q) const = 0;
	//! Get probability distribution function values on BOTTOM boundary
	virtual  std::vector<T> getBottomBoundaryValue(int const q) const = 0;
	//! Get probability distribution function values on LEFT boundary
	virtual std::vector<T> getLeftBoundaryValue(int const q) const = 0;
	//! Get probability distribution function values on RIGHT boundary
	virtual std::vector<T> getRightBoundaryValue(int const q) const = 0;

	//! Set probability distribution function values on TOP boundary equal to parameter array
	virtual void setTopBoundaryValue(int const q, std::vector<T> const & row) = 0;
	//! Set probability distribution function values on BOTTOM boundary equal to parameter array
	virtual void setBottomBoundaryValue(int const q, std::vector<T> const & row) = 0;
	//! Set probability distribution function values on LEFT boundary equal to parameter array
	virtual void setLeftBoundaryValue(int const q, std::vector<T> const & coll) = 0;
	//! Set probability distribution function values on RIGHT boundary equal to parameter array
	virtual void setRightBoundaryValue(int const q, std::vector<T> const & coll) = 0;


#pragma region Methods

	////! Fill each of kQ component of probability distribution function except boundaries with value
	//virtual void fillWithoutBoundaries(T const value) = 0;

	////! Fill each of kQ component of probability distribution function boundaries with value
	//virtual void fillBoundaries(T const value) = 0;

	////! Resize each of kQ component of probability distribution function 
	//virtual void resize(unsigned rows, unsigned colls) = 0;

	////! —читает плотность в кажой из €чеек области
	//virtual MacroscopicParam<T> calculateDensity() const = 0;

	////! —читает скорость в кажой из €чеек области
	//virtual MacroscopicParam<T> calculateVelocity(const double mas[kQ], MacroscopicParam<T> const & density) const = 0;

#pragma endregion





};

#endif // !DISTRIBUTION_FUNC_INTERFACE_H

