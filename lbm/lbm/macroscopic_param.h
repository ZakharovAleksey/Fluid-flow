#pragma once

#include"my_matrix.h"
/*!
	Класс для реализации макроскопических параметров, таких как плотность, скорость.

	Реализован на основе класса матрицы, то есть представляет из себя ПОЛЕ плотности,
	и скорости соответственно.
*/
template<typename T>
class MacroscopicParam : public Matrix<T>
{
public:
	MacroscopicParam();
	MacroscopicParam(unsigned rows, unsigned colls);
	virtual ~MacroscopicParam() {}

private:

};

#include"macroscopic_param_impl.h"