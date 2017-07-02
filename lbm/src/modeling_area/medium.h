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
	// Obstacle inside the modeling area
	OBSTACLE			= 7,
};

//! Node of obstacle in fluid flow
struct ObstacleNode
{
	//! Obstacles coordinates
	const int y_;
	const int x_;

	//! Value which will further store distribution function for given q [1 : kQ]
	double value_;

	ObstacleNode(const int y, const int x, double val) : y_(y), x_(x), value_(val) {}
};


#pragma region 2d

//! Modeling area implementation class.
// Consists from a matrix filled with values determing type of current node
class Medium
{
public:
	Medium();
	Medium(int rows, int colls);
	~Medium();

	//! Return true, selected node represents fluid
	bool IsFluid(const int  y, const int x) const;

	//! Resize current Medium with values  !!! DELETE THIS METHOD AND MAKE USING OINTERS LIKE IN 3D
	void Resize(const int  rows, const int colls);

	// Переписать через метод size() реализованный у класса Matrix<>
	std::pair<unsigned int, unsigned int> size() const;

	//! Returns true if medium has additional immersed bodies in modeling area, but common TOP, BOTTOM, LEFT, RIGTH boundaries
	bool IsIncludeObstacles() const;
	//! Returns number of rows
	const int GetRowsNumber() const;
	//! Returns number of columns
	const int GetColumnsNumber() const;



	void AddCircle(const int x0, const int y0, const int radius);
	void AddSquare(const int leftX, const int leftY, const int width, const int height);
	void AddBottomAngle(const int leftX, const int leftY, const int radius);
	void AddTopAngle(const int leftX, const int leftY, const int radius);
	void AddRightAngle(const int leftX, const int leftY, const int radius);
	void AddCircleTopFalf(const int x0, const int y0, const int radius);


	friend std::ostream & operator<<(std::ostream & os, Medium const & medium);
	//! Returns 
	NodeType Get(const int y, const int x) const;

	//! Gets array of obstacles
	std::vector<ObstacleNode> GetObstacles() const;

private:

	//! Fill modeling ara with inital boundary-fluid parameters: all nodes, but boundaries - fluid, while other
	//! nodes are left, right, top, bottom boundaries
	void FillInitialState();

private:

	//! Number of rows in modeling area
	int rows_;
	//! Number of columns in modeling area
	int colls_;

	//! Stored array with obstacles nodes in modeling area
	std::vector<ObstacleNode> obstacles_;

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

