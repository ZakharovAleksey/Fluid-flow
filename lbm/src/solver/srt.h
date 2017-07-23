#pragma once

#include<memory>
#include<windows.h>

#include <fstream> // file streams
#include <sstream> // string streams

#include"solver.h"
#include"..\modeling_area\fluid.h"
#include"..\modeling_area\medium.h"
#include"bc\bc.h"


#include"im_body\immersed_body.h" // Only for Microphone - next to do in separate class

#pragma region 2d

// Пока тут - но перенести

struct Vector2D
{
	double x;
	double y;

	Vector2D() : x(0), y(0) {}
	Vector2D(double x, double y) : x(x), y(y) {}
};

const Vector2D gravity(0.0, 0.0001);

// SRT approach implementation.
//
// Relaxation parameter tau must be bigger then 0.5 to achive good results.
// It is better to choose it near 1.0;
class SRTsolver : iSolver
{
public:
	SRTsolver() : tau_(0.0) {}
	SRTsolver(double const tau, Medium & medium, Fluid & fluid, BCs* bc);
	virtual ~SRTsolver() {}

	virtual void feqCalculate() override;
	virtual void Streaming() override;
	virtual void Collision() override;

	virtual void Solve(int iteration_number) override;
	virtual void Recalculate() override;

	// FOrce
	void ForceCalculation()
	{
		for (int y = 1; y < fluid_->GetRowsNumber(); ++y)
			for (int x = 1; x < fluid_->GetColumnsNumber(); ++x)
			{
				force_.at(0)(y, x) = (1.0 - 0.5 / tau_) * kW[0] * (3.0 * ((-fluid_->vx_(y, x) * (gravity.x) + -fluid_->vy_(y, x) * (gravity.y))));
				force_.at(1)(y, x) = (1.0 - 0.5 / tau_) * kW[1] * (3.0 * ((1.0 - fluid_->vx_(y, x)) * (gravity.x) + (-fluid_->vy_(y, x)) * (gravity.y)) + 9.0 * (fluid_->vx_(y, x)) * (gravity.x));
				force_.at(2)(y, x) = (1.0 - 0.5 / tau_) * kW[3] * (3.0 * ((-fluid_->vx_(y, x)) * ( gravity.x) + (1 - fluid_->vy_(y, x)) * (gravity.y)) + 9.0 * (fluid_->vy_(y, x)) * (gravity.y));
				force_.at(3)(y, x) = (1.0 - 0.5 / tau_) * kW[2] * (3.0* ((-1.0 - fluid_->vx_(y, x)) * (gravity.x) + (-fluid_->vy_(y, x)) * (gravity.y)) + 9.0 * (fluid_->vx_(y, x)) * (gravity.x));
				force_.at(4)(y, x) = (1.0 - 0.5 / tau_) * kW[4] * (3.0 * ((-fluid_->vx_(y, x)) * (gravity.x) + (-1 - fluid_->vy_(y, x)) * (gravity.y)) + 9.0 * (fluid_->vy_(y, x)) * (gravity.y));
				force_.at(5)(y, x) = (1.0 - 0.5 / tau_) * kW[5] * (3.0 * ((1.0 - fluid_->vx_(y, x)) * (gravity.x) + (1 - fluid_->vy_(y, x)) * (gravity.y)) + 9.0 * (fluid_->vx_(y, x) + fluid_->vy_(y, x)) * (gravity.x + gravity.y));
				force_.at(6)(y, x) = (1.0 - 0.5 / tau_) * kW[8] * (3.0 * ((-1.0 - fluid_->vx_(y, x)) * (gravity.x) + (1 - fluid_->vy_(y, x)) * (gravity.y)) + 9.0 * (fluid_->vx_(y, x) - fluid_->vy_(y, x)) * (gravity.x - gravity.y));
				force_.at(7)(y, x) = (1.0 - 0.5 / tau_) * kW[6] * (3.0 * ((-1.0 - fluid_->vx_(y, x)) * (gravity.x) + (-1 - fluid_->vy_(y, x)) * (gravity.y)) + 9.0 * (fluid_->vx_(y, x) + fluid_->vy_(y, x)) * (gravity.x + gravity.y));
				force_.at(8)(y, x) = (1.0 - 0.5 / tau_) * kW[7] * (3.0 * ((1.0 - fluid_->vx_(y, x)) * (gravity.x) + (-1 - fluid_->vy_(y, x)) * (gravity.y)) + 9.0 * (fluid_->vx_(y, x) - fluid_->vy_(y, x)) * (gravity.x - gravity.y));
			}
	
	}

protected:

	//! Creates folder for output data if not existed yet
	void CreateDataFolder(std::string folder_name) const;

protected:
	//! Relaxation parameter
	double tau_;

	std::array<Matrix2D<double>, kQ> force_;

	Medium* medium_;
	Fluid* fluid_;
	BCs* bc_;
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
	void GetNFProfile(const MacroscopicParam3D<double> & physVal, const int chan_numb, const int iter_numb);

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





