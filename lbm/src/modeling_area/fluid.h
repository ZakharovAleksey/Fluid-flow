#pragma once


#include <fstream> // file streams
#include <sstream> // string streams

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

	Fluid(int rows, int colls, const Medium & medium);
	~Fluid();
	// !!! TRY TO REPLASE EVERYWHERE !!!
	std::pair<unsigned, unsigned> size() const;
	//! Returns number of rows
	const int GetRowsNumber() const;
	//! Returns number of columns
	const int GetColumnsNumber() const;


	void Poiseuille_IC(double const dvx);

	//!  Writes current density and velocity fields to *.vtk file
	void WriteFluidToVTK(const std::string path, const int time);

private:

	//! Fills initial field for density based on medium (takes into account immersed bodies represented in medium)
	void FillInitialDensity(const Medium & medium);

private:
	//! Number of rows in fluid area
	int rows_;
	//! Number of columns in fluid area
	int colls_;

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

	//! Returns depth number of fluid domain, or number size along Z-axis
	int GetDepthNumber() const;
	//! Returns rows number of fluid domain, or number size along Y-axis
	int GetRowsNumber() const;
	//! Returns columns number of fluid domain, or number size along X-axis
	int GetColumnsNumber() const;

	//! Applies Poiseuille initial condition to left boundary
	void PoiseuilleIC(double const dvx);

	//! Set 'q'-s component of distribution function with choosen value
	void SetDistributionFuncValue(const int q, double const value);

	// !!! —делать Matrix<> - const
	// Gets 'q'-s component of distribution finction on layer depth 'z'
	Matrix2D<double> GetDistributionFuncLayer(const int z, const int q);
	// Sets 'q'-s component of distribution finction on layer depth 'z' is equal to 'value'
	void SetDistributionFuncLayerValue(const int z, const int q, const int value);

	// Recalculate dencity for each node in fluid domain
	void RecalculateRho();
	// Racalculate all three velocities: vx, vy, vz for each node in fluid domain
	void RecalculateV();
	// Total rho calculation of all fluid domain (For check onlly)
	long double TotalRho();

	void WriteFluidToVTK(const std::string path, int time);

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



