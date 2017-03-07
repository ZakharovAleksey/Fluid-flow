#pragma once

#include"../phys_values/2d/macroscopic_param_2d.h"
#include"../phys_values/2d/distribution_func_2d.h"
#include"medium.h"

class SRTsolver;

/*!
	Stores all parameters for fluid describing, i.e:
		- Fluid density field
		- Fluid velocity field (Two components)
		- Probability distribution function field
		- Equilibrium probability distribution function field
*/
class Fluid
{
	friend class SRTsolver;

public:


	Fluid(unsigned rows, unsigned colls);
	~Fluid();


	std::pair<unsigned, unsigned> size() const;

	void Poiseuille_IC(double const dvx);



private:

	//! Rows count for fluid modeling area
	unsigned rows_;
	//! Columns count for fluid modeling area
	unsigned colls_;

// ”ЅЅЅЅЅЅ–––јјјјј“№№№№№№№№ - это было просто дл€ тестироани€ чтобы передать f_ ка аргумент дл€ BCs
public:
	//! Fluid density field
	MacroscopicParam<double> rho_;
	//! Fluid velocity X-component field
	MacroscopicParam<double> vx_;
	//! Fluid velocity Y-component field
	MacroscopicParam<double> vy_;

	//! Probability distribution function field
	DistributionFunction<double> f_;
	//! Equilibrium probability distribution function field
	DistributionFunction<double> feq_;

};
