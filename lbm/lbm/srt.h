#pragma once

#include<memory>

#include"solver.h"
#include"fluid.h"
#include"medium.h"

#include"bc.h"


/*!
	Класс реализующий алгоритмы для расчета жидкости в рамках SRT подхода.

	Время релаксации \tau >= 0.5 для сходимости. Лучше выбирать около 1.0
*/
class SRTsolver : iSolver
{
public:
	SRTsolver(double const tau, Medium & medium, Fluid & fluid);
	virtual ~SRTsolver() {}

	virtual void feq_calculate();

	virtual void streaming();
	virtual void collision();

	virtual void solve(int iteration_number);
	virtual void recalculate();

private:
	//! Время релаксации
	double const tau_;

	Medium* medium_;
	Fluid* fluid_;
};

