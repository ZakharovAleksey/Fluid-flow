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

	Fluid(unsigned rows, unsigned colls);
	~Fluid();

	std::pair<unsigned, unsigned> size() const;
	void Poiseuille_IC(double const dvx);

	void AddImmersedBodies(const Medium & medium)
	{
		assert(medium.size().first == rows_);
		assert(medium.size().second == colls_);

		for(int y = 1; y < rows_ - 1; ++y)
			for (int x = 1; x < colls_ - 1; ++x)
			{
				if(medium.Get(y,x) == NodeType::BODY_IN_FLUID)
				{
					rho_(y, x) = 0.0;
					vx_(y, x) = 0.0;
					vy_(y, x) = 0.0;
				}
			}
	}


	void write_fluid_vtk(std::string path,  int time) {

		// Create filename

		path += "/fluid_t" + std::to_string(time) + ".vtk";

		//std::stringstream output_filename;
		//output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
		std::ofstream output_file;

		output_file.open(path); //output_filename.str().c_str());
		if (output_file.is_open())
		{

			// Write VTK header

			int Nx = size().second;
			int Ny = size().first;

			output_file << "# vtk DataFile Version 3.0\n";
			output_file << "fluid_state\n";
			output_file << "ASCII\n";
			output_file << "DATASET RECTILINEAR_GRID\n";
			output_file << "DIMENSIONS " << Nx << " " << Ny - 2 << " 1" << "\n";
			output_file << "X_COORDINATES " << Nx << " float\n";

			for (int X = 0; X < Nx; ++X)
				output_file << X + 0.5 << " ";

			output_file << "\n";
			output_file << "Y_COORDINATES " << Ny - 2 << " float\n";

			for (int Y = 1; Y < Ny - 1; ++Y)
				output_file << Y - 0.5 << " ";

			output_file << "\n";
			output_file << "Z_COORDINATES " << 1 << " float\n";
			output_file << 0 << "\n";
			output_file << "POINT_DATA " << Nx * (Ny - 2) << "\n";

			output_file << "SCALARS density float 1\n";
			output_file << "LOOKUP_TABLE default\n";

			for (int Y = 1; Y < Ny - 1; ++Y)
				for (int X = 0; X < Nx; ++X)
					output_file << rho_(Y, X) << "\n";

			output_file << "SCALARS vx float 1\n";
			output_file << "LOOKUP_TABLE default\n";

			for (int Y = 1; Y < Ny - 1; ++Y)
				for (int X = 0; X < Nx; ++X)
					output_file << vx_(Y, X) << "\n";

			output_file << "SCALARS vy float 1\n";
			output_file << "LOOKUP_TABLE default\n";

			for (int Y = 1; Y < Ny - 1; ++Y)
				for (int X = 0; X < Nx; ++X)
					output_file << vy_(Y, X) << "\n";

			output_file << "VECTORS velocity float\n";

			for (int Y = 1; Y < Ny - 1; ++Y)
				for (int X = 0; X < Nx; ++X)
					output_file << vx_(Y, X) << " " << vy_(Y, X) << " 0\n";
		}
		else
		{
			std::cout << "Could not open file for fluid vtk writing!\n";
		}
		
		output_file.close();
	}

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



