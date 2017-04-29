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

	// Removes current version of the matrix. Allocates memory for new matrix and fills it with zeros.
	virtual void Resize(int new_rows_numb, int new_colls_numb, int new_depth_numb = 0) = 0;
};
