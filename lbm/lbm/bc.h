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
	TOP,
	BOTTOM,
	RIGHT,
	LEFT,
};

typedef std::pair<Boundary, BCType> Wall_info;

class BCs
{
	friend class Fluid;
public:
	BCs(unsigned rows, unsigned colls, DistributionFunction<double> & dfunc);
	~BCs();
	
	/*!
		Записывает необходимые для применения boundary_condition_type ГУ функции распределения
		в соответствующие BC поля класса.

		Пример: get_values(TOP, PERIODIC) - заполнит top_boundary_ 2, 5, 6 компонентой с верхней
		границы.
	*/ 
	bool get_values(Boundary const BC, BCType const boundary_condition_type);
	
	//! Подгатавливает все значения для дальнейшего рассчета ГУ
	//! Порядок заполнения типа ОБЯЗАТЕЛЬНО (TOR, BOTTOM, LEFT, RIGHT)
	void prepair_bc_values(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc);

	//! Заполняет соответствующие компоненты функции распределения вычисленными в BC значениями
	void set_values(Boundary const BC, BCType const boundary_condition_type);
	
	//! 
	void recording_bc_values(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc);

#pragma region different_BC_implementation
	//! Периодические ГУ
	void periodic_bc(Boundary const first, Boundary const second);
	//! ГУ типа отскока
	void bounce_back_bc(Boundary const first);
	//! ГУ типа Фон-Неймана (постоянный поток вдоль поверхности)
	void von_neuman_bc(Boundary const first, Fluid & fluid, double const v);

	void apply_bc(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc);

#pragma endregion

	friend std::ostream & operator<<(std::ostream & os, BCs const & BC);

private:
	//! Меняет в std::map<> first аргумент с from на to не изменяя second
	void swap_id(std::map<int, std::vector<double> > & map, int const from, int const to);

private:
	//! Число строк в матрице, хранящий значения на верхней и нижней стенке [= rows_ матрицы]
	unsigned length_;
	//! Число строк в матрице, хранящей значения на правой и левой стенке [= colls_ - 2 матрицы]
	unsigned height_;

	//! Указатель на исползуемую в Fluid функцию распределения, чтобы работать с границами
	DistributionFunction<double>* f_ptr_;

	//! Хранит индекс компоненты функции распределения и ее значения на верхней стенке
	//!	Пример: top_boundary_[1] = { Значения f[1] на верхней границе }
	std::map<int, std::vector<double> > top_boundary_;
	//! Хранит индекс компоненты функции распределения и ее значения на нижней стенке
	std::map<int, std::vector<double> > bottom_boundary_;
	//! Хранит индекс компоненты функции распределения и ее значения на левой стенке
	std::map<int, std::vector<double> > left_boundary_;
	//! Хранит индекс компоненты функции распределения и ее значения на правой стенке
	std::map<int, std::vector<double> > right_boundary_;
};
