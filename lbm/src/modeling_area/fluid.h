#pragma once

#include"../phys_values/2d/macroscopic_param_2d.h"
#include"../phys_values/2d/distribution_func_2d.h"
#include"../phys_values/3d/macroscopic_param_3d.h"
#include"../phys_values/3d/distribution_func_3d.h"

#include"medium.h"
#include"../solver/solver.h"

class SRTsolver;

#pragma region 2d

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

#pragma endregion

#pragma region 3d

#include<memory>

/*!
Stores all parameters for fluid describing, i.e:
- Fluid density field
- Fluid velocity field (Two components)
- Probability distribution function field
- Equilibrium probability distribution function field
*/


class Fluid3D
{
	friend class SRTsolver;

	typedef std::unique_ptr<MacroscopicParam3D<double>> MacroscopicParamPtr;
	typedef std::unique_ptr<DistributionFunction3D<double>> DistributionFuncPtr;

public:
	Fluid3D(int depth, int rows, int colls);
	~Fluid3D() {}

	int GetDepthNumber() const;
	int GetRowsNumber() const;
	int GetColumnsNumber() const;

	void Poiseuille_IC(double const dvx);

	// Set q-s component of distribution function with choosen value
	void SetDistributionFuncValue(const int q, double const value);

	// !!! —делать Matrix<> - const
	// Get q-s component of distribution finction on layer depth 'z'
	Matrix2D<double> GetDistributionFuncLayer(const int z, const int q);
	// Set q-s component of distribution finction on layer depth 'z' is equal to 'value'
	void SetDistributionFuncLayerValue(const int z, const int q, const int value);

	// Recalculate dencity
	void RecalculateRho();
	// Racalculate all three velocities: vx, vy, vz
	void RecalculateV();
	

private:

	//! Number of rows (Y-axis size  value)
	int rows_;
	//! Number of columns (X-axis size  value)
	int colls_;
	//! Depth (Z-axis size  value)
	int depth_;

	// Recalculate single velocity component
	void RecalculateVelocityComponent(const MacroscopicParamPtr & v_ptr, const int e[]);

public:

	//! Fluid density field
	MacroscopicParamPtr rho_;
	//! Fluid velocity X-component field
	MacroscopicParamPtr vx_;
	//! Fluid velocity Y-component field
	MacroscopicParamPtr vy_;
	//! Fluid velocity Z-component field
	MacroscopicParamPtr vz_;

	//! Probability distribution function field
	DistributionFuncPtr f_;
	//! Equilibrium probability distribution function field
	DistributionFuncPtr feq_;

};

#pragma endregion



