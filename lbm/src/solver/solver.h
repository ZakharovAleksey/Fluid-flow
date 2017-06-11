#pragma once

#include"../phys_values/2d/distribution_func_2d.h"
#include"../phys_values/3d/distribution_func_3d.h"

#pragma region 2d



// Model of velocity directions in D2Q9 model
//  6 2 5  ^ y
//   \|/   |
//  3 0 1  O---> x
//   /|\
//  7 4 8


//! Weigth for probability distribution function calculation
const double kW[kQ]{ 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };

//! X-components witch determ particle movement
const double kEx[kQ]{ 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0 };

//! Y-components witch determ particle movement
const double kEy[kQ]{ 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0 };

//! Y-компоненты определе€ющие направлени€ распространени€ псевдочастиц (умножеди на -1 чтобы up = 0, boottom = rows)
//const double kEy[kQ]{ 0.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0 };

#pragma endregion

#pragma region 3d

#include<vector>

//! This is the vector of weights : initilize in feqCalculation() because it is too big
static void FillWeightsFor3D(std::vector<double> & w)
{
	w.resize(kQ3d, 1.0 / 36.0);

	w.at(0) = 12.0 / 36.0;
	w.at(9) = 2.0 / 36.0;
	w.at(14) = 2.0 / 36.0;

	for (int i = 1; i <= 4; ++i)
		w.at(i) = 2.0 / 36.0;
	
}


//! Directions in assordace with Dmitry Biculov article
const int ex[kQ3d] { 0, 1,   0, -1,   0,  1,   -1, -1,   1, 0,   1,  0,   -1, 0,    0,  1,   0, -1,  0 };
const int ey[kQ3d] { 0, 0,  -1,  0,   1, -1,   -1,  1,   1, 0,   0, -1,    0, 1,    0,  0,  -1,  0,  1 };
const int ez[kQ3d] { 0, 0,   0,  0,   0,  0,    0,  0,   0,-1,  -1, -1,   -1,-1,    1,  1,   1,  1,  1 };

#pragma endregion


//! Interface for different LBM approaches, such as:
// - SRT (single-relaxation-time)
// - MRT (multiple-relaxation-time)
// - IB-SRT (single-relaxation-time with possibility of immersed bodies modeling)
class iSolver
{
public:

	virtual ~iSolver() {}

	//! Performs equilibrium probability distribution function calculation
	virtual void feqCalculate() = 0;

	//! Performs streaming of particles to neighbour nodes
	virtual void Streaming() = 0;
	//! Performs collision of particles in nodes of gread
	virtual void Collision() = 0;

	//! Performs step-by-step modeling of fluid flow
	virtual void Solve(int iteration_number) = 0;

	//! Performs calculation of macroscopic parameters, such as fluid density, velocity
	virtual void Recalculate() = 0;

private:

};
