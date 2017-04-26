#pragma once

#include<memory>


#include"solver.h"
#include"..\modeling_area\fluid.h"
#include"..\modeling_area\medium.h"
#include"bc\bc.h"



#pragma region 2d

/*!
SRT approach implementation.

Relaxation parameter tau must be bigger then 0.5 to achive good results.
It is better to choose it near 1.0;

*/
class SRTsolver : iSolver
{
public:
	SRTsolver(double const tau, Medium & medium, Fluid & fluid);
	virtual ~SRTsolver() {}

	//! Equilibrium probability distribution function calculation in SRT implementation
	virtual void feqCalculate();

	//! Streaming of particles to neighbour nodes in SRT implementation
	virtual void streaming();
	//! Collision of particles in nodes in SRT implementation
	virtual void collision();

	//! Solver for modeling procedure in SRT implementation
	virtual void solve(int iteration_number);

	//! Recalculation procedure (recalculate density, velocity) in SRT implementation
	virtual void recalculate();

private:
	//! Relaxation parameter
	double const tau_;

	Medium* medium_;
	Fluid* fluid_;
};


#pragma endregion


#pragma region 3d

class SRT3DSolver : iSolver
{
public:
	SRT3DSolver(double const tau, Medium3D & medium, Fluid3D & fluid);
	virtual ~SRT3DSolver() {}

	virtual void feqCalculate();

	//! Streaming of particles to neighbour nodes in SRT implementation
	virtual void streaming();
	//! Collision of particles in nodes in SRT implementation
	virtual void collision() {}

	virtual void solve(int iteration_number)
	{
		fluid_->Poiseuille_IC(0.01);

		feqCalculate();
		std::cout << *fluid_->feq_;
	}

	//! Recalculation procedure (recalculate density, velocity) in SRT implementation
	virtual void recalculate() {}

private:
	//! Relaxation parameter
	double const tau_;

	// make them smart pointer

	Medium3D* medium_;
	Fluid3D* fluid_;
};


#pragma endregion



