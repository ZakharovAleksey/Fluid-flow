#pragma once

#ifndef MRT_H
#define MRT_H

#include"srt.h"

#pragma region 2d 

class MRTSolver : SRTsolver
{
public:

	MRTSolver(double const tau, Medium & medium, Fluid & fluid);
	virtual ~MRTSolver() {}

	void Collision() override;
	void Solve(int iteration_number) override;

private:

	// Matrix, which transforms the distribution function f to the velocity moment m (A.A. Mohammad 2012)
	const double M_[kQ][kQ] =
	{
		{ 1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0 },
		{ -4.0,	-1.0,	-1.0,	-1.0,	-1.0,	2.0,	2.0,	2.0,	2.0 },
		{ 4.0,	-2.0,	-2.0,	-2.0,	-2.0,	1.0,	1.0,	1.0,	1.0 },
		{ 0.0,	1.0,	0.0,	-1.0,	0.0,	1.0,	-1.0,	-1.0,	1.0 },
		{ 0.0,	-2.0,	0.0,	2.0,	0.0,	1.0,	-1.0,	-1.0,	1.0 },
		{ 0.0,	0.0,	1.0,	0.0,	-1.0,	1.0,	1.0,	-1.0,	-1.0 },
		{ 0.0,	0.0,	-2.0,	0.0,	2.0,	1.0,	1.0,	-1.0,	-1.0 },
		{ 0.0,	1.0,	-1.0,	1.0,	-1.0,	0.0,	0.0,	0.0,	0.0 },
		{ 0.0,	0.0,	0.0,	0.0,	0.0,	1.0,	-1.0,	1.0,	-1.0 }
	};

	// Inverse matrix, which transforms the velocity moment m, back to the distribution function f (A.A. Mohammad 2012)
	const double Minv_[kQ][kQ] =
	{
		{ 4.0 / 36.0,    -4.0 / 36.0,     4.0 / 36.0,     0.0 / 36.0,     0.0 / 36.0,     0.0 / 36.0,     0.0 / 36.0,     0.0 / 36.0,     0.0 / 36.0, },
		{ 4.0 / 36.0,    -1.0 / 36.0,    -2.0 / 36.0,     6.0 / 36.0,    -6.0 / 36.0,     0.0 / 36.0,     0.0 / 36.0,     9.0 / 36.0,     0.0 / 36.0, },
		{ 4.0 / 36.0,    -1.0 / 36.0,    -2.0 / 36.0,     0.0 / 36.0,     0.0 / 36.0,     6.0 / 36.0,    -6.0 / 36.0,    -9.0 / 36.0,     0.0 / 36.0, },
		{ 4.0 / 36.0,    -1.0 / 36.0,    -2.0 / 36.0,    -6.0 / 36.0,     6.0 / 36.0,     0.0 / 36.0,     0.0 / 36.0,     9.0 / 36.0,     0.0 / 36.0, },
		{ 4.0 / 36.0,    -1.0 / 36.0,    -2.0 / 36.0,     0.0 / 36.0,     0.0 / 36.0,    -6.0 / 36.0,     6.0 / 36.0,    -9.0 / 36.0,     0.0 / 36.0, },
		{ 4.0 / 36.0,     2.0 / 36.0,     1.0 / 36.0,     6.0 / 36.0,     3.0 / 36.0,     6.0 / 36.0,     3.0 / 36.0,     0.0 / 36.0,     9.0 / 36.0, },
		{ 4.0 / 36.0,     2.0 / 36.0,     1.0 / 36.0,    -6.0 / 36.0,    -3.0 / 36.0,     6.0 / 36.0,     3.0 / 36.0,     0.0 / 36.0,    -9.0 / 36.0, },
		{ 4.0 / 36.0,     2.0 / 36.0,     1.0 / 36.0,    -6.0 / 36.0,    -3.0 / 36.0,    -6.0 / 36.0,    -3.0 / 36.0,     0.0 / 36.0,     9.0 / 36.0, },
		{ 4.0 / 36.0,     2.0 / 36.0,     1.0 / 36.0,     6.0 / 36.0,     3.0 / 36.0,    -6.0 / 36.0,    -3.0 / 36.0,     0.0 / 36.0,    -9.0 / 36.0, }
	};

	// Transformation matrix S = diag(1.0, 1.2, 1.0, 1.0, 1.2, 1.0, 1.2, 1/tau, 1/tau) (A.A. Mohammad 2012)
	std::vector<double> S_;

	// Matrix, which is the result of multiplication of S transformation matrix with M^{-1} : MinvS_ = M^{-1} * S
	Matrix2D<double> MinvS_;
};


#pragma endregion

#endif // !MRT_H


