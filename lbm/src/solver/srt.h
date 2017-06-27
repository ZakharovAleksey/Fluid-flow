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

protected:

	//! Creates folder for output data if not existed yet
	void CreateDataFolder(std::string folder_name) const;

protected:
	//! Relaxation parameter
	double tau_;

	Medium* medium_;
	Fluid* fluid_;
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



