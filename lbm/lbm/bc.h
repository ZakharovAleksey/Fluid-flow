#pragma once

#include"fluid.h"
#include<memory>

#include<map>


//! Указатель на DistriputionFunction
typedef std::unique_ptr<DistributionFunction<double>> distr_func_ptr;

enum class BCType {
	PERIODIC,
	BOUNCE_BACK,
	VON_NEUMAN,
};

enum class Boundary {
	UP,
	BOTTOM,
	RIGHT,
	LEFT,
};

class BCs
{
	friend class Fluid;
public:
	BCs(unsigned rows, unsigned colls, DistributionFunction<double> & dfunc);
	~BCs();

	//! Метод, заполняющий поля класса BCs соответствующими значениями функций распределения
	bool get_boundady_values();

	/*!
		Периодические граничные условия.

		Простейший тип ГУ, при котором каждая из kQ омпонент скорости при достижении границы области
	моделирования перемещается на противоположную границу.

	Пример: при достижениии правой грницы компоненты f[1], f[5], f[8] переходят в те-же компоненты
	на левой границе и f[3], f[6], f[7] левой границы переходят в те-же компоненты правой
	границы.

	*/
	void set_periodic_bc(Boundary first, Boundary second);

	/*!
		Граничные условия типа отскока.

		ГУ при которых каждая из kQ компонент функции распределения при достижении границы области
		модлирования заменяется на противоположную.

		Пример: На правой границе f[1] переходит в f[3], f[5]->f[6], f[8]->f[7] на той-же границе.
	*/
	void set_bounce_back_bc(Boundary first);

	/*
		Граничные условия Фон-Неймана.

		ГУ при которых сохраняется постоянным значение потока на выбранной поверхности. Скорость
		потока v передается в качстве аргумента функции. 

		Требует передачу Fluid по ссылке, так как при этом типе ГУ изменяется плотность и скорость.
		НА ДАННЫЙ МОМЕНТ ГУ РЕАЛИЗОВАННО ТОЛЬКО ДЛЯ ЛЕВОЙ СТЕНКИ.
	*/
	void set_von_neumann_bc(Boundary first, Fluid & fluid, double const v);

	friend std::ostream & operator<<(std::ostream & os, BCs const & BC);
private:
	//! Метод, заполняющий одно из полей класса BCs соответствующими значениями функции распределения 
	bool get_values(Boundary BC);

private:
	//! Число строк в матрице, хранящий значения на верхней и нижней стенке [= rows_ матрицы]
	unsigned length_;
	//! Число строк в матрице, хранящей значения на правой и левой стенке [= colls_ - 2 матрицы]
	unsigned height_;

	//! Указатель на исползуемую в Fluid функцию распределения, чтобы работать с границами
	distr_func_ptr f_ptr_;

	/*!
		Массив матриц у которой в q-ой строке лежит q-ая функция распределения на верхней стенке

		f[0]( на верхней границе )
		f[1]( на верхней границе )
		...
		f[8]( на верхней границе )

	*/
	std::array<Matrix<double>, kQ> up_boundary_;

	std::map<unsigned, Matrix<double> > top_;
	

	//!	Массив матриц у которой в q-ой строке лежит q-ая функция распределения на нижней стенке	
	std::array<Matrix<double>, kQ> bottom_boundary_;
	
	//!	Массив матриц у которой в q-ой строке лежит q-ая функция распределения на правой стенке
	std::array<Matrix<double>, kQ> right_boundary_;

	//!	Массив матриц у которой в q-ой строке лежит q-ая функция распределения на левой стенке	
	std::array<Matrix<double>, kQ> left_boundary_;

};
