#pragma once

#include<memory>
#include<map>

#include"../../modeling_area/fluid.h"

//! SmartPoiner to DistriputionFunction
typedef std::unique_ptr<DistributionFunction<double>> distr_func_ptr;

//! Stores Boundary Conditions type index
enum class BCType {
	PERIODIC,
	BOUNCE_BACK,
	VON_NEUMAN,
};

//! Stores Boundary type index
enum class Boundary {
	TOP,
	BOTTOM,
	RIGHT,
	LEFT,
	// In 3d case only
	CLOSE_IN,
	FAAR
};


/*!
	Main idea of BC appling process:
	1. Store all necessary probability distribution function values on chosen boundary before STREAMING.
	2. Change this values depending on choosen BC type : Periodic, Bounce Back, e.t.c.
	3. Record this values to an appropriate probability distribution functions values
*/

#pragma region 2d

class BCs
{
	friend class Fluid;

public:

	BCs(DistributionFunction<double> & dfunc);

	~BCs();
	
	//! Prepare ALL probability distribution function values for BC applying
	//! Порядок заполнения типа ОБЯЗАТЕЛЬНО (TOR, BOTTOM, LEFT, RIGHT)
	void PrepareValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc);
	//! Record ALL boundaries with already calculated BC distribution function values 
	void RecordValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc);


	//! Applies periodic boundary conditions
	void PeriodicBC(Boundary const first, Boundary const second);
	//! Applies bounce back boundary conditions
	void BounceBackBC(Boundary const first);
	//! Applies Von-Neumann boundary conditions (constant velocity flow on choosen boundary)
	void VonNeumannBC(Boundary const first, Fluid & fluid, double const vx, double const vy);

	friend std::ostream & operator<<(std::ostream & os, BCs const & BC);


	void swap(double& a, double & b)
	{
		double temp = a;
		a = b;
		b = a;
	}

	void PrepareAdditionalBCs(Medium & medium)
	{
		for (int q = 0; q < kQ; ++q)
		{
			std::vector<ImmersedBodyVal> v;

			for (int y = 1; y < f_ptr_->size().first - 1; ++y)
			{
				for (int x = 1; x < f_ptr_->size().second - 1; ++x)
				{
					if (medium.Get(y, x) == NodeType::BODY_IN_FLUID && f_ptr_->Get(q, y, x) != 0.0)
					{
						v.push_back(ImmersedBodyVal(y, x, f_ptr_->Get(q, y, x)));
						f_ptr_->Set(q, y, x, 0.0);
					}
				}
			}

			additionalBCs.insert(std::make_pair(q, v));
		}

		/*for (auto i : additionalBCs)
		{
			std::cout << i.first << " : ";
			for (auto j : i.second)
				std::cout << j.distrFuncValue_ << "(" << j.y_ << " , " << j.x_ << ") ";
			std::cout << std::endl;
		}*/
	}

	void AdditionalBounceBackBCs()
	{
		std::swap(additionalBCs.at(1), additionalBCs.at(3));
		std::swap(additionalBCs.at(2), additionalBCs.at(4));
		std::swap(additionalBCs.at(5), additionalBCs.at(7));
		std::swap(additionalBCs.at(6), additionalBCs.at(8));

		/*for (auto i : additionalBCs)
		{
			std::cout << i.first << " : ";
			for (auto j : i.second)
				std::cout << j.distrFuncValue_ << "(" << j.y_ << " , " << j.x_ << ") ";
			std::cout << std::endl;
		}*/
	}

	void RecordAdditionalBCs()
	{
		for (auto row : additionalBCs)
		{
			for (auto j : row.second)
				f_ptr_->Set(row.first, j.y_ + ey[row.first], j.x_ + ex[row.first], j.distrFuncValue_);
		}
		additionalBCs.clear();
	}


private:

	//! Prepare values for CHOOSEN ONE BC BEFORE Streaming
	bool PrepareValuesForSingleBC(Boundary const BC, BCType const boundary_condition_type);
	//! Writes a single distribution function components to an appropriate class field
	bool WriteBoundaryValues(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids, std::vector<double>(DistributionFunction<double>::*ptrToFunc)(int)const);
	bool WriteVonNeumannBoundaryValues(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids_1, const std::vector<int> & bc_ids_2, std::vector<double>(DistributionFunction<double>::*ptrToFunc)(int)const);


	//! Record values for choosen ONE BC AFTER Streaming (BC applying itself)
	void RecordValuesOnSingleBC(Boundary const BC, BCType const boundary_condition_type);
	//! Records a single component distribution function component boundary to appropriate boundary of distribution function
	bool RecordBoundaryValues(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids, void(DistributionFunction<double>::* ptrToFunc)(int const, std::vector<double> const &));

	//! Swap two stored boundary values
	void SwapId(std::map<int, std::vector<double> > & map, int const from, int const to);

	
	//! Directly calculate all distribution function components for choosen boundary
	void CalculateVonNeumanBCValues(Boundary const first, const int size, std::map<int, std::vector<double>> & boundary,
		/* Ids of disr. functions are necessary for calculations: example for top ew need {0,1,3} and {2,5,6} */const std::vector<int> ids_1, const std::vector<int> ids_2,
		double const vx, double const vy);

private:

	//! Ids of approppriate boundaries
	const std::vector<int> top_ids_{ 2,5,6 };
	const std::vector<int> bottom_ids_{ 4,7,8 };
	const std::vector<int> left_ids_{ 3,6,7 };
	const std::vector<int> right_ids_{ 1,5,8 };

	//! Ids of boundaries for Von-Neumann BCs
	const std::vector<int> mid_height_ids_{ 0,2,4 };
	const std::vector<int> mid_width_ids_{ 0,1,3 };

	//! Poiner to Fluid distribution function to work with it's boundaries (Попробовать переделать через ссылку)
	DistributionFunction<double>* f_ptr_;

	//! Store index of probability distribution function and it's values on TOP boundary
	//!	Example: top_boundary_[1] = { Store values f[1] on top boundary }
	std::map<int, std::vector<double> > top_boundary_;
	std::map<int, std::vector<double> > bottom_boundary_;
	std::map<int, std::vector<double> > left_boundary_;	
	std::map<int, std::vector<double> > right_boundary_;


	// for additional BCs

	struct ImmersedBodyVal
	{
		
		int y_;
		int x_;

		double distrFuncValue_;

		ImmersedBodyVal(int y, int x, double value) : y_(y), x_(x), distrFuncValue_(value) {}
	};

	typedef std::map<int, std::vector<ImmersedBodyVal> > VectorOfMap;

	VectorOfMap additionalBCs;

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

	void VonNeumannBC(Boundary const first, const double vx = 0.0, const double vy = 0.0, const double vz = 0.0);

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
	bool WriteVonNeumannBoundaryValues(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids_1, const std::vector<int> & bc_ids_2, std::vector<double>(DistributionFunction3D<double>::*ptrToFunc)(int)const);


	//! Record values for choosen ONE BC AFTER Streaming (BC applying itself)
	bool RecordValuesForSingleBC(Boundary const BC, BCType const boundary_condition_type);
	//! Records a single component distribution function component boundary to appropriate boundary of distribution function
	bool RecordBoundaryValues(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids, void(DistributionFunction3D<double>::* ptrToFunc)(int const, std::vector<double> const &));

	//! Swap two stored boundary values
	void SwapIds(std::map<int, std::vector<double> > & map, int const from, int const to);

private:

	//! Indexes of velocity components on appropriate boundaries
	const std::vector<int> top_ids_{ 9,10,11,12,13 };
	const std::vector<int> bottom_ids_{ 14,15,16,17,18 };
	const std::vector<int> left_ids_{ 3,6,7,12,17 };
	const std::vector<int> right_ids_{ 1,5,8,10,15 };
	const std::vector<int> near_ids_{ 4,7,8,13,18 };
	const std::vector<int> far_ids_{ 2,5,6,11,16 };
	
	// For Von Neumann BC
	const std::vector<int> middle_layer_ids_{ 0,1,2,3,4,5,6,7,8 };

	// !!! Не нужны - их убрать !!!

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