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
	DIRICHLET,
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

#pragma region BC Info

//! Contains all information about boundary : boundary wall, type of bpundary condition, e.t.c.
//! Parent class to all boundary consdition information classes : periodic, bounce-back, Von-Neumann, Dirichlet
class BCInfo
{
public:
	BCInfo(Boundary boundary, BCType bc_type) : wall_(boundary), type_(bc_type) {}

	//! Returns type of boundary condition on chosen wall
	BCType GetType() const { return type_; }
	//! Returns x-component velocity value for Von-Neumann boundary condition, 0 otherwise
	virtual std::vector<double> GetVelocityX() const { return std::vector<double>(); };
	//! Returns y-component velocity value for Von-Neumann boundary condition, 0 otherwise
	virtual std::vector<double> GetVelocityY() const { return std::vector<double>();};
	//! Returns density value for Von-Neumann boundary condition, 0 otherwise
	virtual std::vector<double> GetDensity() const { return std::vector<double>(); };
protected:
	//! Wall of modeling area (TOP, BOTTOM, LEFT, RIGTH)
	Boundary wall_;
	//! Type of boundary condition on this wall
	BCType type_;
};

//! Contains all information about boundary for periodic boundary condition
class BCInfoPeriodic : public BCInfo
{
public:
	BCInfoPeriodic(Boundary first) : BCInfo(first, BCType::PERIODIC) {}
	~BCInfoPeriodic() {}
private:
};

//! Contains all information about boundary for bounce-back boundary condition
class BCInfoBounceBack : public BCInfo
{
public:
	BCInfoBounceBack(Boundary boundary) : BCInfo(boundary, BCType::BOUNCE_BACK) {}
	~BCInfoBounceBack() {}
};

//! Contains all information about boundary for Von-Neumann boundary condition
class BCInfoVonNeumann : public BCInfo
{
public:
	BCInfoVonNeumann(Boundary boundary, std::vector<double> vx, std::vector<double> vy) : BCInfo(boundary, BCType::VON_NEUMAN), vx_(vx), vy_(vy) {}
	~BCInfoVonNeumann() {}

	std::vector<double> GetVelocityX() const override { return vx_; }
	std::vector<double> GetVelocityY() const override { return vy_; }
private:
	//! Velocities along x and y axis respectively for chosen wall
	std::vector<double> vx_;
	std::vector<double> vy_;
};

//! Contains all information about boundary for Dirichlet boundary condition
class BcInfoDirichlet : public BCInfo
{
public:
	BcInfoDirichlet(Boundary boundary, std::vector<double> rho) : BCInfo(boundary, BCType::DIRICHLET), rho_(rho) {}
	~BcInfoDirichlet() {}

	std::vector<double> GetDensity() const override { return rho_; }
private:
	//! Dencsity for chosen boundary wall
	std::vector<double> rho_;
};

#pragma endregion

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

	//! Set periodic boundary conditions to chosen boundaries
	void SetPeriodicBC(Boundary first, Boundary second);
	//! Set bounce-back boundary conditions to chosen boundary
	void SetBounceBackBC(Boundary boundary);
	//! Set Von-Neumann boundary conditions to chosen boundary
	void SetVonNeumannBC(Boundary boundary, std::vector<double> vx, std::vector<double> vy);
	//! Set Dirichlet boundary conditions to chosen boundary
	void SetDirichletBC(Boundary boundary, std::vector<double> rho);

	//! Prepares values of distribution function for all boundaries and obstacles, if last takes place
	void PrepareBCValues(const Medium & medium);
	//! Records values of distribution function for all boundaries and obstacles, if last takes place
	void RecordBCValues(const Medium & medium);

	//! Performs chosen earlier boundary conditions 
	void PerformBC(Fluid* fluid, const Medium & medium);

	friend std::ostream & operator<<(std::ostream & os, BCs const & BC);

private:

	//! Prepare ALL probability distribution function values for BC applying
	//! Порядок заполнения типа ОБЯЗАТЕЛЬНО (TOR, BOTTOM, LEFT, RIGHT)
	void PrepareValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc, const Medium & medium);
	//! Record ALL boundaries with already calculated BC distribution function values 
	void RecordValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc, const Medium & medium);

	//! Prepares values for bounce-back boundary conditions for obstacles in fluid
	void PrepareValuesForObstacles(const Medium & medium);
	//! Performs bounce-back boundary conditions for obstacles in fluid
	void ObstaclesBounceBackBCs();
	//! Records values for Bounce-back boundary conditions for obstacles in fluid
	void RecordValuesForObstacles();

	//! Applies periodic boundary conditions
	void PeriodicBC(Boundary const first, Boundary const second);
	//! Applies bounce back boundary conditions
	void BounceBackBC(Boundary const first);
	//! Applies Von-Neumann boundary conditions
	void VonNeumannBC(Boundary const first, Fluid & fluid, std::vector<double> const vx, std::vector<double> const vy);
	//! Applies Dirichlet boundary conditions
	void DirichletBC(Boundary const first, Fluid & fluid, std::vector<double> const rho_0);

	//! Prepare values for CHOOSEN ONE BC BEFORE Streaming
	bool PrepareValuesForSingleBC(Boundary const BC, BCType const boundary_condition_type);
	//! Writes a single distribution function components to an appropriate class field (BC: Periodic, Bounce-Back )
	bool WriteBoundaryValuesFirstOrder(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids, std::vector<double>(DistributionFunction<double>::*ptrToFunc)(int)const);
	//! Writes a single distribution function components to an appropriate class field (BC: Von-Neumann, Dirichle )
	bool WriteBoundaryValuesSecondOrder(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids_1, const std::vector<int> & bc_ids_2, std::vector<double>(DistributionFunction<double>::*ptrToFunc)(int)const);


	//! Record values for choosen ONE BC AFTER Streaming (BC applying itself)
	void RecordValuesOnSingleBC(Boundary const BC, BCType const boundary_condition_type);
	//! Records a single component distribution function component boundary to appropriate boundary of distribution function
	bool RecordBoundaryValues(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids, void(DistributionFunction<double>::* ptrToFunc)(int const, std::vector<double> const &));

	//! Swap two stored boundary values
	void SwapId(std::map<int, std::vector<double> > & map, int const from, int const to);

	
	//! Directly calculate all distribution function components for choosen boundary
	void CalculateVonNeumanBCValues(Boundary const first, const int size, std::map<int, std::vector<double>> & boundary,
		/* Ids of disr. functions are necessary for calculations: example for top ew need {0,1,3} and {2,5,6} */const std::vector<int> ids_1, const std::vector<int> ids_2,
		std::vector<double> const vx, std::vector<double> const vy);

	//! Directly calculate all distribution function components for choosen boundary
	void CalculateDirichletBCValues(Boundary const first, const int size, std::map<int, std::vector<double>> & boundary,
		/* Ids of disr. functions are necessary for calculations: example for top ew need {0,1,3} and {2,5,6} */const std::vector<int> ids_1, const std::vector<int> ids_2,
		std::vector<double> const rho_0);

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

	//! Stores index of probability distribution function and it's values on TOP boundary
	//!	Example: top_boundary_[1] = { Store values f[1] on top boundary }
	std::map<int, std::vector<double> > top_boundary_;
	std::map<int, std::vector<double> > bottom_boundary_;
	std::map<int, std::vector<double> > left_boundary_;	
	std::map<int, std::vector<double> > right_boundary_;

	//! Stores index of probability distribution function and it's values on obstacles nodes
	std::map<int, std::vector<ObstacleNode>> obstacles_bc_;

	//! Performs boundary conditions (Bounce-back, Von-Neumann or Dirichlet) which could be applied to 
	//! one (single) boundary
	void PerformSingleBoundaryBC(std::shared_ptr<BCInfo> bc_info, Boundary boundary, Fluid* fluid)
	{
		switch (bc_info->GetType())
		{
		case BCType::PERIODIC:
			// In this case do nothing because periodic boundary conditions applies ONLY to TWO boundaries
			// and algorithm will deal with them later
			break;
		case BCType::BOUNCE_BACK:
			BounceBackBC(boundary);
			break;
		case BCType::VON_NEUMAN:
			VonNeumannBC(boundary, *fluid, bc_info->GetVelocityX(), bc_info->GetVelocityY());
			break;
		case BCType::DIRICHLET:
			DirichletBC(boundary, *fluid, bc_info->GetDensity());
			break;
		default:
			std::cout << "Error! You are tryeing to perform wrong boundary conditions to ONE boundary!\n";
			break;
		}
	}

	//! Pointer which stores information for boundary conditions on each boundary
	std::shared_ptr<BCInfo> top_info_;
	std::shared_ptr<BCInfo> bottom_info_;
	std::shared_ptr<BCInfo> right_info_;
	std::shared_ptr<BCInfo> left_info_;
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
	const std::vector<int> top_ids_{ 4,8,12,14,17 };//{ 9,10,11,12,13 };
	const std::vector<int> bottom_ids_{5,9,13,15,16};//{ 14,15,16,17,18 };
	const std::vector<int> left_ids_{ 1,10,11,12,13 };//{ 3,6,7,12,17 };
	const std::vector<int> right_ids_{0,6,7,8,9 };//{ 1,5,8,10,15 };
	const std::vector<int> near_ids_{3,7,11,16,17 };//{ 4,7,8,13,18 };
	const std::vector<int> far_ids_{ 2,6,10,14,15 }; //{ 2,5,6,11,16 };
	
	// For Von Neumann BC
	const std::vector<int> middle_layer_ids_{ 0, 1, 2, 3, 6, 7, 10, 11, 18 };

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