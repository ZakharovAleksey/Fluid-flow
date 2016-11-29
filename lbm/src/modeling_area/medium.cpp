#include"medium.h"
/*!
	Конструктор класса medium:

	Матрица состоящая из NodeType.
	По умолчанию область моделирования заполняется жидкостью : NodeType::FLUID;
		- после этого по умолчанию задается верхняя, нижняя, левая, правая границы соответственно 
		NodeType::UPPER_BOUNDARY, NodeType::BOTTOM_BOUNDARY, NodeType::LEFT_BOUNDARY, NodeType::RIGHT_BOUNDARY.

	- в дальнейшем необходимо проверять область на наличие погруженных тел.

	Пример:
	1	1	1	1
	3	0	0	4
	3	0	0	4
	2	2	2	2

*/

Medium::Medium() : rows_(0), colls_(0), medium_() {}


Medium::Medium(int rows, int colls):
	rows_(rows), colls_(colls) 
{ 

	assert(rows_ > 2 && colls_ > 2);

	medium_.Resize(rows_, colls_);

	for (int x = 0; x < colls_; ++x) 
	{
		medium_(0, x) = NodeType::UPPER_BOUNDARY;
		medium_(rows_ - 1, x) = NodeType::BOTTOM_BOUNDARY;
	}

	for (int y = 1; y < rows_ - 1; ++y) 
	{
		medium_(y, 0) = NodeType::LEFT_BOUNDARY;
		medium_(y, colls_ - 1) = NodeType::RIGHT_BOUNDARY;
	}
}

Medium::~Medium() {}

bool Medium::is_fluid(int y, int x) const
{
	assert(y < rows_ && x < colls_);
	return medium_(y, x) == NodeType::FLUID;
}

void Medium::resize(int rows, int colls)
{
	rows_ = rows;
	colls_ = colls;

	medium_.Resize(rows_, colls_);

	// Потом выделить это в отдельную функцию и реализовать через нее конструктор и resize()
	assert(rows_ > 2 && colls_ > 2);

	medium_.Resize(rows_, colls_);

	for (int x = 0; x < colls_; ++x) 
	{
		medium_(0, x) = NodeType::UPPER_BOUNDARY;
		medium_(rows_ - 1, x) = NodeType::BOTTOM_BOUNDARY;
	}

	for (int y = 1; y < colls_ - 2; ++y) 
	{
		medium_(y, 0) = NodeType::LEFT_BOUNDARY;
		medium_(y, colls_ - 1) = NodeType::RIGHT_BOUNDARY;
	}
}


std::pair<int, int> Medium::size() const
{
	return std::make_pair(rows_, colls_);
}


std::ostream & operator<<(std::ostream & os, Medium const & medium) {

	for (int y = 0; y < medium.rows_; ++y) {
		for (int x = 0; x < medium.colls_; ++x) {

			if (medium.medium_(y, x) == NodeType::FLUID)
				os << std::setw(3) << 0;
			else if (medium.medium_(y, x) == NodeType::UPPER_BOUNDARY)
				os << std::setw(3) << 1;
			else if (medium.medium_(y, x) == NodeType::BOTTOM_BOUNDARY)
				os << std::setw(3) << 2;
			else if (medium.medium_(y, x) == NodeType::LEFT_BOUNDARY)
				os << std::setw(3) << 3;
			else if (medium.medium_(y, x) == NodeType::RIGHT_BOUNDARY)
				os << std::setw(3) << 4;
		}
		os << std::endl;
	}

	return os;
}