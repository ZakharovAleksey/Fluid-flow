#pragma once

#include<memory>

#include"solver.h"
#include"fluid.h"
#include"medium.h"

#include"bc.h"

typedef std::unique_ptr<Medium> Medium_ptr;
typedef std::unique_ptr<Fluid> Fluid_ptr;

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

	virtual void solve();

private:
	//! Время релаксации
	double const tau_;

	Medium_ptr medium_;
	Fluid_ptr fluid_;
};

