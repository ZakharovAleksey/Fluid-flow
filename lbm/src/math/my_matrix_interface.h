#pragma once

#include<vector>

// Interface for Matrix class.
template<typename T>
class iMatrix
{
public:
	
	virtual ~iMatrix() {}

	// Returns the sum of all elements of the matrix
	virtual long double GetSum() const = 0;
	
	// Returns vector, with values from row with choosen index from matrix
	virtual std::vector<T> GetRow(unsigned const y) const = 0;

	// Set elements of "y" row in the current matrix equal to values from std::vector "row"
	virtual void SetRow(unsigned const y, std::vector<T> const & row) = 0;
	
	// Returns std::vector, with values in range [1 : rows_ - 2] from column with "x" index.
	virtual std::vector<T> GetColumn(unsigned const x) const = 0;

	// Set elements of "x" column in the current matrix equal to values from std::vector "coll"
	virtual void SetColumn(unsigned const x, std::vector<T> const & coll) = 0;

private:

};
