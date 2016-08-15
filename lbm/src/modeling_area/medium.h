#pragma once

#include <type_traits>

#include"..\math\my_matrix.h"

/*!
	Тип ячейки: жидкость, или ячейка на границе области моделирования.
	В дальнейшем расширение для типа погруженной границы.

	Реализован как enum class чтобы избежать возможные конфликты имен.
*/
enum class NodeType : int 
{
	FLUID				= 0,
	UPPER_BOUNDARY		= 1,
	BOTTOM_BOUNDARY		= 2,
	LEFT_BOUNDARY		= 3,
	RIGHT_BOUNDARY		= 4,
};

/*!
	Класс реализующий область моделирования.

	Представляет собой матрицу NodeType-ов в которой отмечена информация о типе ячейки.
*/
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

	void resize(unsigned rows, unsigned colls);

	friend std::ostream & operator<<(std::ostream & os, Medium const & medium);

	// Переписать через метод size() реализованный у класса Matrix<>
	std::pair<unsigned int, unsigned int> size() const;

private:
	unsigned rows_;
	unsigned colls_;

	Matrix<NodeType> medium_;
};
