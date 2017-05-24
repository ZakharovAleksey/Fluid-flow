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


Medium::Medium() : rows_(0), colls_(0), medium_() {}


Medium::Medium(unsigned rows, unsigned colls) :
	rows_(rows), colls_(colls)
{
	assert(rows_ > 2 && colls_ > 2);

	medium_.Resize(rows_, colls_);

	for (int x = 0; x < colls_; ++x)
	{
		medium_(0, x) = NodeType::TOP_BOUNDARY;
		medium_(rows_ - 1, x) = NodeType::BOTTOM_BOUNDARY;
	}

	for (int y = 1; y < rows_ - 1; ++y)
	{
		medium_(y, 0) = NodeType::LEFT_BOUNDARY;
		medium_(y, colls_ - 1) = NodeType::RIGHT_BOUNDARY;
	}
}

Medium::~Medium() {}

bool Medium::is_fluid(unsigned y, unsigned x) const
{
	assert(y < rows_ && x < colls_);
	return medium_(y, x) == NodeType::FLUID;
}

void Medium::resize(unsigned rows, unsigned colls)
{
	rows_ = rows;
	colls_ = colls;

	medium_.Resize(rows_, colls_);

	// Потом выделить это в отдельную функцию и реализовать через нее конструктор и resize()
	assert(rows_ > 2 && colls_ > 2);

	medium_.Resize(rows_, colls_);

	for (int x = 0; x < colls_; ++x)
	{
		medium_(0, x) = NodeType::TOP_BOUNDARY;
		medium_(rows_ - 1, x) = NodeType::BOTTOM_BOUNDARY;
	}

	for (int y = 1; y < colls_ - 2; ++y)
	{
		medium_(y, 0) = NodeType::LEFT_BOUNDARY;
		medium_(y, colls_ - 1) = NodeType::RIGHT_BOUNDARY;
	}
}

std::pair<unsigned, unsigned> Medium::size() const
{
	return std::make_pair(rows_, colls_);
}

void Medium::AddCircleInMedium(const int x0, const int y0, const int radius)
{
	assert(x0 + radius < colls_ - 1 && x0 - radius > 1);
	assert(y0 + radius < rows_ - 1 && y0 - radius > 1);

	int xStart = x0 - radius; int xStop = x0 + radius;
	int yStart = y0 - radius; int yStop = y0 + radius;

	for (int y = yStart; y < yStop; ++y)
	{
		for (int x = xStart; x < xStop; ++x)
		{
			if ( pow(x- x0, 2) + pow(y - y0, 2) < pow(radius, 2))
				medium_(y, x) = NodeType::BODY_IN_FLUID;
		}
	}
}

void Medium::AddCircleTopFalf(const int x0, const int y0, const int radius)
{
	assert(x0 + radius < colls_ - 1 && x0 - radius > 1);
	assert(y0 < rows_ - 1 && y0 >= 1);

	int xStart = x0 - radius; int xStop = x0 + radius;
	int yStart = y0; int yStop = y0 + radius;

	for (int y = yStart; y < yStop; ++y)
	{
		for (int x = xStart; x < xStop; ++x)
		{
			if (pow(x - x0, 2) + pow(y - y0, 2) < pow(radius, 2))
				medium_(y, x) = NodeType::BODY_IN_FLUID;
		}
	}
}

void Medium::AddCircleBottomFalf(const int x0, const int y0, const int radius)
{
	assert(x0 + radius < colls_ - 1 && x0 - radius > 1);
	assert(y0 <= rows_ - 1 && y0 >= 1);

	int xStart = x0 - radius; int xStop = x0 + radius;
	int yStart = y0 - radius; int yStop = y0;

	for (int y = yStart; y < yStop; ++y)
	{
		for (int x = xStart; x < xStop; ++x)
		{
			if (pow(x - x0, 2) + pow(y - y0, 2) < pow(radius, 2))
				medium_(y, x) = NodeType::BODY_IN_FLUID;
		}
	}
}


std::ostream & operator<<(std::ostream & os, Medium const & medium) {

	for (int y = 0; y < medium.rows_; ++y) {
		for (int x = 0; x < medium.colls_; ++x) {

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
			else if (medium.medium_(y, x) == NodeType::BODY_IN_FLUID)
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