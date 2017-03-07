#pragma once

#include"../phys_values/2d/distribution_func_2d.h"

//! Weigth for probability distribution function calculation
const double kW[kQ]{ 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
 				   1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };


/*
	Directions in D2Q9 model
	  6   2   5
	   \  |  /
	3 --  0  -- 1
	   /  |  \
	  7   4   8
*/

//! X-components witch determ particle movement
const double kEx[kQ]{ 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0 };

//! Y-components witch determ particle movement
const double kEy[kQ]{ 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0 };

//! Y-компоненты определе€ющие направлени€ распространени€ псевдочастиц (умножеди на -1 чтобы up = 0, boottom = rows)
//const double kEy[kQ]{ 0.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0 };

/*!
	Interface for different simulation aproaches:
		- SRT
		- MRT
		- Immersed Boundary SRT
*/
class iSolver
{
public:

	virtual ~iSolver() {}

#pragma region Methods

	//! Equilibrium probability distribution function calculation
	virtual void feqCalculate() = 0;

	//! Streaming of particles to neighbour nodes
	virtual void streaming() = 0;
	//! Collision of particles in nodes
	virtual void collision() = 0;

	//! Solver for modeling procedure
	virtual void solve(int iteration_number) = 0;

	//! Recalculation procedure (recalculate density, velocity)
	virtual void recalculate() = 0;

#pragma endregion

private:

};
