#pragma once

#include<vector>

/// <summary>
/// Interface for Matrix class.
/// </summary>
/// <remarks>
/// This class provides all necessary function for future LBM methods implementation depending
/// on the type of executed task - 2 dimensional task or 3-dimensional.
/// </remarks>
template<typename T>
class iMatrix
{
public:
	
	virtual ~iMatrix() {}

	/// <summary>
	/// Returns the sum of all elements of the matrix
	/// </summary>
	/// <returns></returns>
	virtual long double GetSum() const = 0;

	/// <summary>
	/// Returns vector, with values from row with choosen index from matrix
	/// </summary>
	/// <param name="y"> Index of row </param>
	/// <returns> vector with values from choosen row </returns>
	virtual std::vector<T> GetRow(unsigned const y) const = 0;

	/// <summary>
	/// Set elements of "y" row in the current matrix equal to values from std::vector "row"
	/// </summary>
	/// <param name="y"> Index of matrix row for insertion operation. </param>
	/// <param name="row"> Vector for for insertion operation. </param>
	virtual void SetRow(unsigned const y, std::vector<T> const & row) = 0;

	/// <summary>
	/// Returns std::vector, with values in range [1 : rows_ - 2] from column with "x" index.
	/// </summary>
	/// <param name="x"> Index of column from witch method extract values. </param>
	/// <returns> Vector filled with values from apropriate column of the current matrix. </returns>
	virtual std::vector<T> GetColumn(unsigned const x) const = 0;

	/// <summary>
	/// Set elements of "x" column in the current matrix equal to values from std::vector "coll"
	/// </summary>
	/// <param name="x"> Index of matrix column for insertion operation. </param>
	/// <param name="coll">  Vector for for insertion operation. </param>
	virtual void SetColumn(unsigned const x, std::vector<T> const & coll) = 0;

private:

};
