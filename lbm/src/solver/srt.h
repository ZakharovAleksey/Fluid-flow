#pragma once

#include<memory>
#include<windows.h>

#include <fstream> // file streams
#include <sstream> // string streams

#include"solver.h"
#include"..\modeling_area\fluid.h"
#include"..\modeling_area\medium.h"
#include"bc\bc.h"

#pragma region 2d


// SRT approach implementation.
//
// Relaxation parameter tau must be bigger then 0.5 to achive good results.
// It is better to choose it near 1.0;
class SRTsolver : iSolver
{
public:
	SRTsolver() : tau_(0.0) {}
	SRTsolver(double const tau, Medium & medium, Fluid & fluid);
	virtual ~SRTsolver() {}

	virtual void feqCalculate() override;
	virtual void Streaming() override;
	virtual void Collision() override;

	virtual void Solve(int iteration_number) override;
	virtual void Recalculate() override;

private:

	//! Creates folder for output data if not existed yet
	void CreateDataFolder(std::string folder_name) const;

protected:
	//! Relaxation parameter
	double tau_;

	Medium* medium_;
	Fluid* fluid_;
};


#pragma endregion

#pragma region 2d mrt

class MRTSolver : SRTsolver
{
public:

	MRTSolver(double const tau, Medium & medium, Fluid & fluid) : SRTsolver(tau, medium, fluid)
	{
		std::cout << " --- Input parameters :\n";
		std::cout << "nu = " << (tau - 0.5) / 3.0 << std::endl;
		std::cout << "Re = " << 0.01 * medium_->size().second / ((tau - 0.5) / 3.0) << std::endl;

		// Create folders for next data uploading
		CreateDataFolder("Data");
		CreateDataFolder("Data\\mrt_lbm_data");
		CreateDataFolder("Data\\mrt_lbm_data\\2d");
		CreateDataFolder("Data\\mrt_lbm_data\\2d\\fluid_txt");
		CreateDataFolder("Data\\mrt_lbm_data\\2d\\fluid_vtk");

		// Fill transformation matrix S
		S_ = { 1.0, 1.2, 1.0, 1.0, 1.2, 1.0, 1.2, 1.0 / tau_, 1.0 / tau_ };

		// Allocate memory for S * M^{-1} matrix
		MinvS_.Resize(kQ, kQ);

		// S - diagonal matrix
		Matrix2D<double> S(kQ, kQ);
		for (int i = 0; i < kQ; ++i)
			S(i, i) = S_.at(i);

		// Perform multiplication: S * M^{-1}
		for(int i = 0; i < kQ; ++i)
			for (int j = 0; j < kQ; ++j)
			{
				MinvS_(i, j) = 0.0;
				for (int k = 0; k < kQ; ++k)
					MinvS_(i, j) += Minv_[i][k] * S(k,j);
			}

	}

	virtual ~MRTSolver() {}

	void Collision() override 
	{
		// Obtain domain size
		const int y_size = medium_->size().first;
		const int x_size = medium_->size().second;

		// Distribution function in momentum space
		DistributionFunction<double> dm(y_size, x_size);
		for (int q = 0; q < kQ; ++q)
			dm[q].FillWith(0.0);

		// Calculate m = M * f
		for (int k = 0; k < kQ; ++k)
			for (int m = 0; m < kQ; ++m)
				dm[k] += M_[k][m] * fluid_->f_[m];

		// Performs calculations of values necessary for meq calculation
		Matrix2D<double> rvx = fluid_->rho_.ScalarMultiplication(fluid_->vx_);
		Matrix2D<double> rvy = fluid_->rho_.ScalarMultiplication(fluid_->vy_);
		Matrix2D<double> vxSq = fluid_->vx_.ScalarMultiplication(fluid_->vx_);
		Matrix2D<double> vySq = fluid_->vy_.ScalarMultiplication(fluid_->vy_);
		Matrix2D<double> vSq = vxSq + vySq;

		// Performs dm = m - m_eq : Additional check m_eq = M * f_eq
		dm[0] -= fluid_->rho_;
		dm[1] -= fluid_->rho_.ScalarMultiplication((-2.0 + 3.0 * vSq));
		dm[2] -= fluid_->rho_.ScalarMultiplication((- 3.0 * vSq + 1));
		dm[3] -= rvx;
		dm[4] += rvx;
		dm[5] -= rvy;
		dm[6] += rvy;
		dm[7] -= fluid_->rho_.ScalarMultiplication(vxSq - vySq);
		dm[8] -= fluid_->rho_.ScalarMultiplication(fluid_->vx_.ScalarMultiplication(fluid_->vy_));

		// Performs f(x + vdt, t + dt) = f(x, t) - M^{-1}S * dm
		for (int k = 0; k < kQ; ++k)
			for (int m = 0; m < kQ; ++m)
				fluid_->f_[k] -= MinvS_(k, m) * dm[m];

	}

	void Solve(int iteration_number) override 
	{
		feqCalculate();
		for (int q = 0; q < kQ; ++q)
			fluid_->f_[q] = fluid_->feq_[q];

		BCs BC(fluid_->f_);

		for (int iter = 0; iter < iteration_number; ++iter)
		{
			Collision();
			BC.PrepareValuesForAllBC(BCType::VON_NEUMAN, BCType::VON_NEUMAN, BCType::VON_NEUMAN, BCType::VON_NEUMAN);

			Streaming();

			BC.VonNeumannBC(Boundary::TOP, *fluid_, 0.01, 0.0);
			BC.VonNeumannBC(Boundary::BOTTOM, *fluid_, -0.01, 0.0);
			BC.VonNeumannBC(Boundary::LEFT, *fluid_, 0.0, -0.01);
			BC.VonNeumannBC(Boundary::RIGHT, *fluid_, 0.0, +0.01);

			BC.RecordValuesForAllBC(BCType::VON_NEUMAN, BCType::VON_NEUMAN, BCType::VON_NEUMAN, BCType::VON_NEUMAN);

			Recalculate();

			feqCalculate();

			std::cout << iter << " Total rho = " << fluid_->rho_.GetSum() << std::endl;

			if (iter % 25 == 0)
			{
				Matrix2D<double> v = CalculateModulus(fluid_->vx_, fluid_->vy_);
				v.WriteFieldToTxt("Data\\mrt_lbm_data\\2d\\fluid_txt", "v", iter);
				//fluid_->vx_.WriteFieldToTxt("Data\\mrt_lbm_data\\2d\\fluid_txt", "vx", iter);
				fluid_->write_fluid_vtk("Data\\mrt_lbm_data\\2d\\fluid_vtk", iter);
			}

		}

	}

private:

	//! Creates folder for output data if not existed yet
	void CreateDataFolder(std::string folder_name) const
	{
		// Get path to current directory
		char buffer[MAX_PATH];
		GetModuleFileName(NULL, buffer, MAX_PATH);

		std::string::size_type pos = std::string(buffer).find_last_of("\\/");
		std::string path = std::string(buffer).substr(0, pos);
		path = path.substr(0, path.size() - 6) + "\\" + folder_name;

		char *cstr = new char[path.length() + 1];
		strcpy(cstr, path.c_str());

		// Create folder if not exist yet
		if (GetFileAttributes(cstr) == INVALID_FILE_ATTRIBUTES)
			CreateDirectory(cstr, NULL);
	}

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



#pragma region 3d

// SRT approach implementation in 3D case.
//
// Relaxation parameter tau must be bigger then 0.5 to achive good results.
// It is better to choose it near 1.0;
class SRT3DSolver : iSolver
{
public:

	SRT3DSolver(double tau, Medium3D & medium, Fluid3D & fluid);
	virtual ~SRT3DSolver() {}

	void feqCalculate() override;
	void Streaming() override;
	void Collision() override;
	void Recalculate() override;
	void Solve(int iteration_number) override;


	void GetProfile(const int chan_numb, const int iter_numb);
	//! Implements correct hetmap writing in file ('length' is a number of elements in one line)
	bool WriteHeatMapInFile(const std::string & file_name, const std::vector<double> & data, const int lenght);


private:
	//! Implementation of streaming for 0-8 velocity directions
	void SubStreamingMiddle(const int depth, const int rows, const int colls);
	//! Implementation of streaming for 9-13 velocity directions
	void SubStreamingTop(const int depth, const int rows, const int colls);
	//! Implementation of streaming for 14-18 velocity directions
	void SubStreamingBottom(const int depth, const int rows, const int colls);

	//! Creates folder for output data if not existed yet
	void CreateDataFolder(std::string folder_name) const;

private:
	//! Relaxation parameter
	double const tau_;

	// !!! Make them smart pointer
	//! Medium domain of simulation
	Medium3D* medium_;
	//! Fluid domain of simulation
	Fluid3D* fluid_;
};


#pragma endregion



