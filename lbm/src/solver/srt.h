#pragma once

#include<memory>


#include"solver.h"
#include"..\modeling_area\fluid.h"
#include"..\modeling_area\medium.h"
#include"bc\bc.h"



/*!
	SRT approach implementation.

	Relaxation parameter tau must be bigger then 0.5 to achive good results. 
	It is better to choose it near 1.0;

*/
class SRTsolver : iSolver
{
public:

#pragma region Constructor

	SRTsolver(double const tau, Medium & medium, Fluid & fluid);
	virtual ~SRTsolver() {}

#pragma endregion

#pragma region Methods
	
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

#pragma endregion

private:

#pragma region Fields

	//! Relaxation parameter
	double const tau_;

	Medium* medium_;
	Fluid* fluid_;

#pragma endregion

};

