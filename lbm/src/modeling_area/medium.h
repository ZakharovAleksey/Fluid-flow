#pragma once

#include <type_traits>

#include"..\math\my_matrix_2d.h"

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

/// <summary>
/// Класс реализующий область моделирования.

/// Представляет собой матрицу NodeType - ов в которой отмечена информация о типе ячейки.
/// </summary>
class Medium
{
public:
	Medium();
	Medium(int rows, int colls);
	~Medium();
	
	bool is_fluid(int y, int x) const;
	// Возможно понадобится определить
	/*bool is_upper_boundary(int y, int x) const;
	bool is_bottom_boundary(int y, int x) const;
	bool is_left_boundary(int y, int x) const;
	bool is_right_boundary(int y, int x) const;*/

	/// <summary>
	/// Resize current Medium with values !!!!!
	/// </summary>
	/// <param name="rows"> lol </param>
	/// <param name="colls"> hah </param>
	void resize(int rows, int colls);

	friend std::ostream & operator<<(std::ostream & os, Medium const & medium);

	// Переписать через метод size() реализованный у класса Matrix<>
	std::pair<int, int> size() const;

private:
	/// <summary>
	/// Rows in modaling area
	/// </summary>
	int rows_;
	/// <summary>
	/// Columns in modaling area
	/// </summary>
	int colls_;

	Matrix2D<NodeType> medium_;
};
