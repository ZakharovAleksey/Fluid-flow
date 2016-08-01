#pragma once

#include"macroscopic_param.h"
#include"distribution_func.h"
#include"medium.h"

class SRTsolver;

/*!
	 ласс который хранит в себе все харрактеристики жидкости:
	 - плотность
	 - 2-е компоненты скорости
	 - функцию распределени€ kQ компонент
	 - равновесную функцию распределени€ kQ компонент
*/
class Fluid
{
	friend class SRTsolver;
public:
	Fluid(unsigned rows, unsigned colls);
	~Fluid();

	void Poiseuille_IC(double const dvx);

	std::pair<unsigned, unsigned> size() const;

private:
	unsigned rows_;
	unsigned colls_;

// ”ЅЅЅЅЅЅ–––јјјјј“№№№№№№№№ - это было просто дл€ тестироани€ чтобы передать f_ ка аргумент дл€ BCs
public:
	MacroscopicParam<double> rho_;
	MacroscopicParam<double> vx_;
	MacroscopicParam<double> vy_;


	DistributionFunction<double> f_;
	DistributionFunction<double> feq_;
};
