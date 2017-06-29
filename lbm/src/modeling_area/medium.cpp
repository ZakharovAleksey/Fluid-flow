#include"medium.h"



/*!
	Medium class constructor (idea is identical to 3D, but with some additions):

	Main body of class, is a matrix, which contains data about all nodes types.
	By default all matrix filfilled with fluid type : NodeType::FLUID;
		- After this procedure complite, we set TOP, BOTTOM, LEFT and RIGHT (NEAR and FAR in 3D case) using
		NodeType::UPPER_BOUNDARY, NodeType::BOTTOM_BOUNDARY, NodeType::LEFT_BOUNDARY, NodeType::RIGHT_BOUNDARY

	2D example :

	1	1	1	1
	3	0	0	4
	3	0	0	4
	2	2	2	2

*/

#pragma region 2d


Medium::Medium() : rows_(0), colls_(0), is_immersed_bodies_(false), medium_() {}


Medium::Medium(int rows, int colls) : rows_(rows), colls_(colls), is_immersed_bodies_(false)
{
	assert(rows_ > 2 && colls > 2);

	medium_.Resize(rows_, colls_);

	FillInitialState();
}

Medium::~Medium() {}

bool Medium::IsFluid(const int  y, const int x) const
{
	assert(y < rows_ && x < colls_);
	return medium_(y, x) == NodeType::FLUID;
}

void Medium::Resize(const int  newRows, const int newColls)
{
	assert(newRows > 2 && newColls > 2);

	rows_ = newRows;
	colls_ = newColls;
	medium_.Resize(rows_, colls_);

	FillInitialState();
}

std::pair<unsigned, unsigned> Medium::size() const
{
	return std::make_pair(rows_, colls_);
}

bool Medium::IsImmersedBodies() const
{
	return is_immersed_bodies_;
}

const int Medium::GetRowsNumber() const
{
	return rows_;
}

const int Medium::GetColumnsNumber() const
{
	return colls_;
}

void Medium::AddCircle(const int x0, const int y0, const int radius)
{
	is_immersed_bodies_ = true;
	//assert(x0 + radius < colls_ - 1 && x0 - radius > 1);
	//assert(y0 + radius < rows_ - 1 && y0 - radius > 1);

	int xStart = x0 - radius; int xStop = x0 + radius;
	int yStart = y0 - radius; int yStop = y0 + radius;

	for (int y = yStart; y < yStop; ++y)
	{
		for (int x = xStart; x < xStop; ++x)
		{
			if (y >= 0 && y < rows_ && x > 0 && x < colls_)
				if ( pow(x- x0, 2) + pow(y - y0, 2) < pow(radius, 2))
					medium_(y, x) = NodeType::OBSTACLE;
		}
	}
}

void Medium::AddSquare(const int leftX, const int leftY, const int width, const int height)
{
	is_immersed_bodies_ = true;
	
	for(int y = leftY; y > leftY - height; --y)
		for (int x = leftX; x < leftX + width; ++x)
		{
			if( y >= 0 && y < rows_ && x > 0 && x < colls_)
				medium_(y, x) = NodeType::OBSTACLE;
		}


}

void Medium::AddBottomAngle(const int leftX, const int leftY, const int radius)
{
	is_immersed_bodies_ = true;
	//assert(leftX > 0 && leftX + radius < colls_);
	//assert(leftY - radius > 0 && leftY < rows_);

	int x0 = leftX + radius;
	int y0 = leftY - radius;
	int circleRad = radius - 1;

	for (int y = leftY; y > leftY - radius; --y)
		for (int x = leftX; x < leftX + radius; ++x)
		{
			if (y >= 0 && y < rows_ && x > 0 && x < colls_)
				if ((x - x0) * (x - x0) + (y - y0) * (y - y0) >= circleRad * circleRad)
					medium_(y, x) = NodeType::OBSTACLE;
		}


}

void Medium::AddTopAngle(const int leftX, const int leftY, const int radius)
{
	is_immersed_bodies_ = true;
	//assert(leftX > 0 && leftX + radius < colls_);
	//assert(leftY >= 0 && leftY + radius < rows_);

	int x0 = leftX + radius;
	int y0 = leftY + radius;
	int circleRad = radius - 1;

	for (int y = leftY; y < leftY + radius; ++y)
		for (int x = leftX; x < leftX + radius; ++x)
		{
			if (y >= 0 && y < rows_ && x > 0 && x < colls_)
				if ((x - x0) * (x - x0) + (y - y0) * (y - y0) >= circleRad * circleRad)
					medium_(y, x) = NodeType::OBSTACLE;
		}

}

void Medium::AddCircleTopFalf(const int x0, const int y0, const int radius)
{
	is_immersed_bodies_ = true;
	//assert(x0 + radius < colls_ - 1 && x0 - radius > 1);
	//assert(y0 + radius < rows_ - 1 && y0 - radius > 1);

	int xStart = x0 - radius; int xStop = x0 + radius;
	int yStart = y0 - radius; int yStop = y0;

	for (int y = yStart; y < yStop; ++y)
	{
		for (int x = xStart; x < xStop; ++x)
		{
			if (y >= 0 && y < rows_ && x > 0 && x < colls_)
				if (pow(x - x0, 2) + pow(y - y0, 2) < pow(radius, 2))
					medium_(y, x) = NodeType::OBSTACLE;
		}
	}
}

void Medium::AddCircleBottomFalf(const int x0, const int y0, const int radius)
{
	is_immersed_bodies_ = true;

	assert(x0 + radius < colls_ - 1 && x0 - radius > 1);
	assert(y0 <= rows_ - 1 && y0 >= 1);

	int xStart = x0 - radius; int xStop = x0 + radius;
	int yStart = y0 - radius; int yStop = y0;

	for (int y = yStart; y < yStop; ++y)
	{
		for (int x = xStart; x < xStop; ++x)
		{
			if (pow(x - x0, 2) + pow(y - y0, 2) < pow(radius, 2))
				medium_(y, x) = NodeType::OBSTACLE;
		}
	}
}

void Medium::FillInitialState()
{
	for (int x = 0; x < colls_; ++x)
	{
		medium_(rows_ - 1, x) = NodeType::TOP_BOUNDARY;
		medium_(0, x) = NodeType::BOTTOM_BOUNDARY;
	}

	for (int y = 1; y < rows_ - 1; ++y)
	{
		medium_(y, 0) = NodeType::LEFT_BOUNDARY;
		medium_(y, colls_ - 1) = NodeType::RIGHT_BOUNDARY;
	}
}


std::ostream & operator<<(std::ostream & os, Medium const & medium) {

	for (int y = medium.rows_ - 1; y >= 0; --y) 
	{
		for (int x = 0; x < medium.colls_; ++x) 
		{

			if (medium.medium_(y, x) == NodeType::FLUID)
				os << std::setw(3) << 0;
			else if (medium.medium_(y, x) == NodeType::TOP_BOUNDARY)
				os << std::setw(3) << 1;
			else if (medium.medium_(y, x) == NodeType::BOTTOM_BOUNDARY)
				os << std::setw(3) << 2;
			else if (medium.medium_(y, x) == NodeType::LEFT_BOUNDARY)
				os << std::setw(3) << 3;
			else if (medium.medium_(y, x) == NodeType::RIGHT_BOUNDARY)
				os << std::setw(3) << 4;
			else if (medium.medium_(y, x) == NodeType::OBSTACLE)
				os << std::setw(3) << 7;
		}
		os << std::endl;
	}

	return os;
}

#pragma endregion

#pragma region 3d


Medium3D::Medium3D() : depth_(0), rows_(0), colls_(0), medium_(nullptr) { }

Medium3D::Medium3D(int depth, int rows, int colls) : depth_(depth), rows_(rows), colls_(colls)
{
	medium_ = std::make_unique<Matrix3D<NodeType>>(depth_, rows_, colls_);
	FillMedium();
}

bool Medium3D::IsFluid(int z, int y, int x) const
{
	return (medium_->operator()(z, y, x) == NodeType::FLUID) ? true : false;
}

void Medium3D::Resize(int depth, int rows, int colls)
{
	depth_ = depth;
	rows_ = rows;
	colls_ = colls;

	medium_.release();
	medium_ = std::make_unique<Matrix3D<NodeType>>(depth_, rows_, colls_);

	FillMedium();
}

int Medium3D::GetDepthNumber() const
{
	return depth_;
}

int Medium3D::GetRowsNumber() const
{
	return rows_;
}

int Medium3D::GetColumnsNumber() const
{
	return colls_;
}

void Medium3D::FillMedium()
{
	// TOP and BOTTOM layers of modeling cube (Oxy plane)
	for (int y = 0; y < rows_; ++y)
		for (int x = 0; x < colls_; ++x)
		{
			(*medium_)(0, y, x)	= NodeType::BOTTOM_BOUNDARY;
			(*medium_)(depth_ - 1, y, x) = NodeType::TOP_BOUNDARY;
		}

	for (int z = 1; z < depth_ - 1; ++z)
	{
		// LEFT and RIGHT layers of modeling cube (Oxz plane)
		for (int y = 0; y < rows_; ++y)
		{
			(*medium_)(z, y, 0) = NodeType::RIGHT_BOUNDARY;
			(*medium_)(z, y, colls_ - 1) = NodeType::LEFT_BOUNDARY;
		}

		// NEAR and FAR layers of modeling cube (Oyz plane)
		for (int x = 0; x < colls_; ++x)
		{
			(*medium_)(z, 0, x) = NodeType::FAR_BOUNDARY;
			(*medium_)(z, rows_ - 1, x) = NodeType::NEAR_BOUNDARY;
		}
	}
}

std::ostream & operator<<(std::ostream & os, const Medium3D & m)
{
	for (int z = 0; z < m.depth_; ++z)
	{
		std::cout << "depth : " << z << " ----------- \n";

		for (int y = 0; y < m.rows_; ++y) {
			for (int x = 0; x < m.colls_; ++x) {


				if (m.medium_->operator()(z, y, x) == NodeType::FLUID)
					os << std::setw(3) << 0;
				else if (m.medium_->operator()(z, y, x) == NodeType::TOP_BOUNDARY)
					os << std::setw(3) << 1;
				else if (m.medium_->operator()(z, y, x) == NodeType::BOTTOM_BOUNDARY)
					os << std::setw(3) << 2;
				else if (m.medium_->operator()(z, y, x) == NodeType::LEFT_BOUNDARY)
					os << std::setw(3) << 3;
				else if (m.medium_->operator()(z, y, x) == NodeType::RIGHT_BOUNDARY)
					os << std::setw(3) << 4;
				else if (m.medium_->operator()(z, y, x) == NodeType::NEAR_BOUNDARY)
					os << std::setw(3) << 5;
				else if (m.medium_->operator()(z, y, x) == NodeType::FAR_BOUNDARY)
					os << std::setw(3) << 6;
			}
			os << std::endl;
		}
	}

	return os;
}

#pragma endregion