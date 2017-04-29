#pragma once

#include <type_traits>

#include"../math/2d/my_matrix_2d.h"
#include"../math/3d/my_matrix_3d.h"

//! Type of Eulerian grid node: fluid or one kind of boundary
enum class NodeType : int 
{
	FLUID				= 0,
	TOP_BOUNDARY		= 1,
	BOTTOM_BOUNDARY		= 2,
	LEFT_BOUNDARY		= 3,
	RIGHT_BOUNDARY		= 4,
	// In case of 3D LBM additional boundaries
	NEAR_BOUNDARY		= 5,
	FAR_BOUNDARY		= 6,
};

#pragma region 2d

//! Modeling area implementation class.
// Consists from a matrix filled with values determing type of current node
class Medium
{
public:
	Medium();
	Medium(unsigned rows, unsigned colls);
	~Medium();

	bool is_fluid(unsigned y, unsigned x) const;

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

	//! Checks if current node is fluid
	bool IsFluid(int z, int y, int x) const;
	//! Resize current medium body [As far as I remember do not implemented in this code because of using std::unique_ptr<>]
	void Resize(int depth, int rows, int colls);

	friend std::ostream & operator<<(std::ostream & os, Medium3D const & m);

	//! Returns depth number of medium, or number size along Z-axis
	int GetDepthNumber() const;
	//! Returns rows number of medium, or number size along Y-axis
	int GetRowsNumber() const;
	//! Returns depth number of medium, or number size along X-axis
	int GetColumnsNumber() const;

private:

	//! Fills each node of medium body with appropriate boundary type
	void FillMedium();

private:
	//! Depth (Z-axis size  value)
	int depth_;
	//! Number of rows (Y-axis size  value)
	int rows_;
	//! Number of columns (X-axis size  value)
	int colls_;
	
	typedef std::unique_ptr<Matrix3D<NodeType>> Matrix3DPtr;

	//! Smart pointer to 3D Matrix, representing mediums body
	Matrix3DPtr medium_;
};

#pragma endregion

