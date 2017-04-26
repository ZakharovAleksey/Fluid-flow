#pragma once

#include <type_traits>

#include"../math/2d/my_matrix_2d.h"
#include"../math/3d/my_matrix_3d.h"

// Type of Eulerian grid node: fluid or boundary
// В дальнейшем расширение для типа погруженной границы.

enum class NodeType : int 
{
	FLUID				= 0,
	UPPER_BOUNDARY		= 1,
	BOTTOM_BOUNDARY		= 2,
	LEFT_BOUNDARY		= 3,
	RIGHT_BOUNDARY		= 4,
};

#pragma region 2d

// Modeling area implementation class.
// Consists from a matrix filled with values determing type of current node
class Medium
{
public:
	Medium();
	Medium(unsigned rows, unsigned colls);
	~Medium();

	bool is_fluid(unsigned y, unsigned x) const;
	// Возможно понадобится определить
	/*bool is_upper_boundary(unsigned y, unsigned x) const;
	bool is_bottom_boundary(unsigned y, unsigned x) const;
	bool is_left_boundary(unsigned y, unsigned x) const;
	bool is_right_boundary(unsigned y, unsigned x) const;*/

	/// <summary>
	/// Resize current Medium with values !!!!!
	/// </summary>
	/// <param name="rows"> lol </param>
	/// <param name="colls"> hah </param>
	void resize(unsigned rows, unsigned colls);

	friend std::ostream & operator<<(std::ostream & os, Medium const & medium);

	// Переписать через метод size() реализованный у класса Matrix<>
	std::pair<unsigned int, unsigned int> size() const;

private:
	/// <summary>
	/// Rows in modaling area
	/// </summary>
	unsigned rows_;
	/// <summary>
	/// Columns in modaling area
	/// </summary>
	unsigned colls_;

	Matrix2D<NodeType> medium_;
};

#pragma endregion


#pragma region 3d

#include<memory>

// Modeling area implementation class.
// Consists from a matrix filled with values determing type of current node
class Medium3D
{
public:
	Medium3D();
	Medium3D(int depth, int rows, int colls);
	~Medium3D() {}

	bool IsFluid(int z, int y, int x) const;
	void Resize(int depth, int rows, int colls);

	friend std::ostream & operator<<(std::ostream & os, Medium3D const & m);

	int GetDepthNumber() const;
	int GetRowsNumber() const;
	int GetColumnsNumber() const;

private:

	// Fill each node of medium body with correct type
	void FillMedium();

private:
	//! Number of rows (Y-axis size  value)
	int rows_;
	//! Number of columns (X-axis size  value)
	int colls_;
	//! Depth (Z-axis size  value)
	int depth_;

	typedef std::unique_ptr<Matrix3D<NodeType>> Matrix3DPtr;

	Matrix3DPtr medium_;
};

#pragma endregion

