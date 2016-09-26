#pragma once

#include<memory>
#include<map>

#include"../../modeling_area/fluid.h"

//! SmartPoiner to DistriputionFunction
typedef std::unique_ptr<DistributionFunction<double>> distr_func_ptr;

//! Store Boundary Conditions type index
enum class BCType {
	PERIODIC,
	BOUNCE_BACK,
	VON_NEUMAN,
};

//! Store Boundary type index
enum class Boundary {
	TOP,
	BOTTOM,
	RIGHT,
	LEFT,
};

class BCs
{
	friend class Fluid;

public:

#pragma region Constructor

	BCs(unsigned rows, unsigned colls, DistributionFunction<double> & dfunc);
	~BCs();

#pragma endregion

#pragma region Methods

	/*!
		Store all needed probability distribution function values on chosen boundary before BC is applying in
		aproppriate class field.

		Example: 
			prepareValues(TOP, PERIODIC) - fill top_boundary_ with 2, 5, 6 component of probability distribution function from
			TOP boundary
	*/ 
	bool prepareValuesOnCurrentBoundary(Boundary const BC, BCType const boundary_condition_type);
	
	//! Prepare ALL probability distribution function values for BC applying
	//! Порядок заполнения типа ОБЯЗАТЕЛЬНО (TOR, BOTTOM, LEFT, RIGHT)
	void prepareValuesForBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc);
	
	//! Record appropriate boundary with already calculated BC distribution function values 
	void recordValuesOnCurrentBoundary(Boundary const BC, BCType const boundary_condition_type);	
	
	//! Record ALL boundaries with already calculated BC distribution function values 
	void recordValuesForBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc);

#pragma region Boundary Conditions (BCs) implementation

	//! Periodic Boundary conditions
	void periodicBC(Boundary const first, Boundary const second);

	//! Bounce Back Boundary conditions
	void bounceBackBC(Boundary const first);

	//! Von-Neumann Boundary conditions
	//! Пока не возвращем массив плотностей так как с ним меньшие погрешности
	void vonNeumannBC(Boundary const first, Fluid & fluid, double const vx, std::vector<double> & velocity_x);

#pragma endregion

#pragma endregion

	friend std::ostream & operator<<(std::ostream & os, BCs const & BC);

private:
	//! Меняет в std::map<> first аргумент с from на to не изменяя second
	void swap_id(std::map<int, std::vector<double> > & map, int const from, int const to);

private:

#pragma region Fields

	//! Rows length [equal to rows_ of matrix]  beacuse all nodes takes placr in BC
	unsigned length_;
	//! Columns height [equal to colls_ - 2 of matrix] because UP and DOWN nodes are already counted in TOP and BOTTOM BC
	unsigned height_;

	//! Poiner to Fluid distribution function to work with it's boundaries (Попробовать переделать через ссылку)
	DistributionFunction<double>* f_ptr_;

	//! Store index of probability distribution function and it's values on TOP boundary
	//!	Пример: top_boundary_[1] = { Значения f[1] на верхней границе }
	std::map<int, std::vector<double> > top_boundary_;
	//! Store index of probability distribution function and it's values on BOTTOM boundary
	std::map<int, std::vector<double> > bottom_boundary_;
	//! Store index of probability distribution function and it's values on LEFT boundary
	std::map<int, std::vector<double> > left_boundary_;
	//! Store index of probability distribution function and it's values on RIGHT boundary
	std::map<int, std::vector<double> > right_boundary_;

#pragma endregion

};
