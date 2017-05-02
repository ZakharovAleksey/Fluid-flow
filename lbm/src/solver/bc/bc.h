#pragma once

#include<memory>
#include<map>

#include"../../modeling_area/fluid.h"

//! SmartPoiner to DistriputionFunction
typedef std::unique_ptr<DistributionFunction<double>> distr_func_ptr;

//! Store Boundary Conditions type index
enum class BCType {
	PERIODIC,
	BOUNCE_BACK,
	VON_NEUMAN,
};

//! Store Boundary type index
enum class Boundary {
	TOP,
	BOTTOM,
	RIGHT,
	LEFT,
	// In 3d case only
	NEAR,
	FAAR
};


#pragma region 2d



class BCs
{
	friend class Fluid;

public:

	BCs(unsigned rows, unsigned colls, DistributionFunction<double> & dfunc);
	~BCs();

	/*!
		Store all needed probability distribution function values on chosen boundary before BC is applying in
		aproppriate class field.

		Example: 
			prepareValues(TOP, PERIODIC) - fill top_boundary_ with 2, 5, 6 component of probability distribution function from
			TOP boundary
	*/ 
	bool prepareValuesOnCurrentBoundary(Boundary const BC, BCType const boundary_condition_type);
	
	//! Prepare ALL probability distribution function values for BC applying
	//! Порядок заполнения типа ОБЯЗАТЕЛЬНО (TOR, BOTTOM, LEFT, RIGHT)
	void prepareValuesForBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc);
	
	//! Record appropriate boundary with already calculated BC distribution function values 
	void recordValuesOnCurrentBoundary(Boundary const BC, BCType const boundary_condition_type);	
	
	//! Record ALL boundaries with already calculated BC distribution function values 
	void recordValuesForBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc);



	//! Periodic Boundary conditions
	void periodicBC(Boundary const first, Boundary const second);

	//! Bounce Back Boundary conditions
	void bounceBackBC(Boundary const first);

	//! Von-Neumann Boundary conditions
	//! Пока не возвращем массив плотностей так как с ним меньшие погрешности
	void vonNeumannBC(Boundary const first, Fluid & fluid, double const vx, std::vector<double> & velocity_x);

	friend std::ostream & operator<<(std::ostream & os, BCs const & BC);

private:
	//! Меняет в std::map<> first аргумент с from на to не изменяя second
	void swap_id(std::map<int, std::vector<double> > & map, int const from, int const to);

private:

	//! Rows length [equal to rows_ of matrix]  beacuse all nodes takes placr in BC
	unsigned length_;
	//! Columns height [equal to colls_ - 2 of matrix] because UP and DOWN nodes are already counted in TOP and BOTTOM BC
	unsigned height_;

	//! Poiner to Fluid distribution function to work with it's boundaries (Попробовать переделать через ссылку)
	DistributionFunction<double>* f_ptr_;

	//! Store index of probability distribution function and it's values on TOP boundary
	//!	Пример: top_boundary_[1] = { Значения f[1] на верхней границе }
	std::map<int, std::vector<double> > top_boundary_;
	//! Store index of probability distribution function and it's values on BOTTOM boundary
	std::map<int, std::vector<double> > bottom_boundary_;
	//! Store index of probability distribution function and it's values on LEFT boundary
	std::map<int, std::vector<double> > left_boundary_;
	//! Store index of probability distribution function and it's values on RIGHT boundary
	std::map<int, std::vector<double> > right_boundary_;

};

#pragma endregion


#pragma region 3d


// Class of 3D boundary conditions:
// - Calculate and apply choosen boundary condidtions to each boundary of modeling fluid domain.

class BCs3D
{
	friend class Fluid;

public:

	BCs3D(int rows, int colls, DistributionFunction3D<double> & dfunc);
	~BCs3D() {}

	//! Prepare ALL values for choosen BC BEFORE Streaming
	void PrepareValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc, BCType const near_bc, BCType far_bc);
	//! Record ALL values for choosen BC AFTER Streaming (BC applying itself)
	void RecordValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc, BCType const near_bc, BCType far_bc);
	
	
	//! Applies periodic boundary conditions
	void PeriodicBC(Boundary const first, Boundary const second);
	//! Applies bounce back boundary conditions
	void BounceBackBC(Boundary const first);

	friend std::ostream & operator<<(std::ostream & os, BCs3D const & BC)
	{
		os.precision(3);

		os << "TOP BOUNDARY ------ \n";
		for (auto i : BC.top_boundary_) {
			os << "f[" << i.first << "] = ";
			for (auto j : i.second)
				os << j << ' ';
			os << std::endl;
		}

		os << "BOTTOM BOUNDARY ------ \n";
		for (auto i : BC.bottom_boundary_) {
			os << "f[" << i.first << "] = ";
			for (auto j : i.second)
				os << j << ' ';
			os << std::endl;
		}

		os << "RIGHT BOUNDARY ------ \n";
		for (auto i : BC.right_boundary_) {
			os << "f[" << i.first << "] = ";
			for (auto j : i.second)
				os << j << ' ';
			os << std::endl;
		}

		os << "LEFT BOUNDARY ------ \n";
		for (auto i : BC.left_boundary_) {
			os << "f[" << i.first << "] = ";
			for (auto j : i.second)
				os << j << ' ';
			os << std::endl;
		}

		os << "NEAR BOUNDARY ------ \n";
		for (auto i : BC.near_boundary_) {
			os << "f[" << i.first << "] = ";
			for (auto j : i.second)
				os << j << ' ';
			os << std::endl;
		}

		os << "FAR BOUNDARY ------ \n";
		for (auto i : BC.far_boundary_) {
			os << "f[" << i.first << "] = ";
			for (auto j : i.second)
				os << j << ' ';
			os << std::endl;
		}

		return os;
	}

private:

	//! Prepare values for CHOOSEN ONE BC BEFORE Streaming
	bool PrepareValuesForSingleBC(Boundary const BC, BCType const bc_type);
	//! Writes a single distribution function components to an appropriate class field
	bool WriteBoundaryValues(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids, std::vector<double>(DistributionFunction3D<double>::*ptrToFunc)(int)const);

	//! Record values for choosen ONE BC AFTER Streaming (BC applying itself)
	bool RecordValuesForSingleBC(Boundary const BC, BCType const boundary_condition_type);
	//! Records a single component distribution function component boundary to appropriate boundary of distribution function
	bool RecordBoundaryValues(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids, void(DistributionFunction3D<double>::* ptrToFunc)(int const, std::vector<double> const &));

	//! Swap two stored boundary values
	void swap_id(std::map<int, std::vector<double> > & map, int const from, int const to);

private:

	//! Indexes of velocity components on appropriate boundaries
	const std::vector<int> top_ids_{ 9,10,11,12,13 };
	const std::vector<int> bottom_ids_{ 14,15,16,17,18 };
	const std::vector<int> left_ids_{ 3,6,7,12,17 };
	const std::vector<int> right_ids_{ 1,5,8,10,15 };
	const std::vector<int> near_ids_{ 4,7,8,13,18 };
	const std::vector<int> far_ids_{ 2,5,6,11,16 };
	

	//! Columns height [equal to colls_ - 2 of matrix] because UP and DOWN nodes are already counted in TOP and BOTTOM BC
	unsigned height_;
	//! Rows length [equal to rows_ of matrix]  beacuse all nodes takes placr in BC
	unsigned length_;
	

	//! Poiner to Fluid distribution function to work with it's boundaries (Попробовать переделать через ссылку)
	DistributionFunction3D<double>* f_ptr_;

	//! Store index of probability distribution function and it's values on TOP boundary
	//!	Example: top_boundary_[1] = { Store values f[1] on top boundary }
	std::map<int, std::vector<double> > top_boundary_;
	std::map<int, std::vector<double> > bottom_boundary_;
	std::map<int, std::vector<double> > left_boundary_;
	std::map<int, std::vector<double> > right_boundary_;
	std::map<int, std::vector<double> > near_boundary_;
	std::map<int, std::vector<double> > far_boundary_;
};

#pragma endregion