#pragma once

#include<memory>
#include<windows.h>

#include <fstream> // file streams
#include <sstream> // string streams

#include"solver.h"
#include"..\modeling_area\fluid.h"
#include"..\modeling_area\medium.h"
#include"bc\bc.h"

#define 	M_PI   3.14159265358979323846
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code

#pragma region 2d

/*!
SRT approach implementation.

Relaxation parameter tau must be bigger then 0.5 to achive good results.
It is better to choose it near 1.0;

*/
class SRTsolver : iSolver
{
public:
	SRTsolver(double const tau, Medium & medium, Fluid & fluid);
	virtual ~SRTsolver() {}

	//! Equilibrium probability distribution function calculation in SRT implementation
	virtual void feqCalculate();


	//! Streaming of particles to neighbour nodes in SRT implementation
	virtual void Streaming();
	//! Collision of particles in nodes in SRT implementation
	virtual void Collision();

	//! Solver for modeling procedure in SRT implementation
	virtual void Solve(int iteration_number);

	//! Recalculation procedure (recalculate density, velocity) in SRT implementation
	virtual void Recalculate();

private:
	//! Relaxation parameter
	double const tau_;

	Medium* medium_;
	Fluid* fluid_;
};


#pragma endregion


#pragma region ib-lbm

////! Point in 2D space
//struct Point
//{
//	double x_;
//	double y_;
//
//	Point() : x_(0.0), y_(0.0) {}
//	Point(double x, double y) : x_(x), y_(y) {}
//};

////! Node of immersed in fluid body with
//struct IBNode
//{
//	//! Current node position
//	Point cur_pos_;
//	//! Reference node position
//	Point ref_pos_;
//
//	//! Node velocity along x-axis
//	double vx_;
//	//! Node velocity along y-axis
//	double vy_;
//	//! Elastic force, acting on the node along x-axis
//	double Fx_;
//	//! Elastic force, acting on the node along y-axis
//	double Fy_;
//
//	IBNode() : cur_pos_(), ref_pos_(), vx_(0.0), vy_(0.0), Fx_(0.0), Fy_(0.0) {}
//};

//class ImmersedBody
//{
//public:
//
//	ImmersedBody(int domainX, int domainY, int nodesNumber, Point center, double radius) : domain_x_(domainX), domain_y_(domainY), nodes_num(nodesNumber), center_(center), radius_(radius)
//	{
//		body_.resize(nodes_num, IBNode());
//
//		for (int id = 0; id < nodes_num; ++id)
//		{
//			// Figure of immersed body
//			/*body_.at(id).cur_pos_.x_ = center_.x_ + radius_ * sin(2.0 * M_PI * (double)id / nodes_num);
//			body_.at(id).ref_pos_.x_ = body_.at(id).cur_pos_.x_;
//
//			body_.at(id).cur_pos_.y_ = center_.x_ + radius_ * cos(2.0 * M_PI * (double)id / nodes_num);
//			body_.at(id).ref_pos_.y_ = body_.at(id).cur_pos_.x_;*/
//
//			body_.at(id).cur_pos_.y_ = center.y_ + radius * sin(2. * M_PI * (double)id / nodes_num);
//			body_.at(id).ref_pos_.y_ = center.y_ + radius * sin(2. * M_PI * (double)id / nodes_num);
//			body_.at(id).cur_pos_.x_ = radius * cos(2. * M_PI * (double)id / nodes_num);
//
//			// Parametrization of the red blood cell shape in 2D
//
//			if (body_.at(id).cur_pos_.x_ > 0) 
//			{
//				
//				body_.at(id).cur_pos_.x_ = center.x_ + sqrt(1 - SQ( (center.y_ - body_.at(id).cur_pos_.y_) / radius) ) * (0.207 + 2.00 * SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius) - 1.12 * SQ(SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius))) * radius / 2;
//				body_.at(id).ref_pos_.x_ = center.x_ + sqrt(1 - SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius)) * (0.207 + 2.00 * SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius) - 1.12 * SQ(SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius))) * radius / 2;
//			}
//			else {
//				body_.at(id).cur_pos_.x_ = center.x_ - sqrt(1 - SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius)) * (0.207 + 2.00 * SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius) - 1.12 * SQ(SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius))) * radius / 2;
//				body_.at(id).ref_pos_.x_ = center.x_ - sqrt(1 - SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius)) * (0.207 + 2.00 * SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius) - 1.12 * SQ(SQ((center.y_ - body_.at(id).cur_pos_.y_) / radius))) * radius / 2;
//			}
//
//		}
//
//	}
//	~ImmersedBody() {}
//
//
//	void CalculateForces()
//	{
//		// Set initial values of force equal to zero
//		for (auto & node : body_)
//		{
//			node.Fx_ = 0.0;
//			node.Fy_ = 0.0;
//		}
//
//		const double arcLen = 2.0 * M_PI * radius_ / nodes_num;
//
//		CalculateStrainForces();
//		CalculateBendingForces();
//	}
//
//	void SpreadForces(Matrix2D<double> &fx, Matrix2D<double> & fy)
//	{
//		fx.FillWith(0.0);
//		fy.FillWith(0.0);
//
//		for (int i = 0; i < nodes_num; ++i)
//		{
//			int x_int = (int)(body_.at(i).cur_pos_.x_ - 0.5 + domain_x_) - domain_x_;
//			int y_int = (int)(body_.at(i).cur_pos_.y_ + 0.5);
//
//			for (int X = x_int; X <= x_int + 1; ++X) {
//				for (int Y = y_int; Y <= y_int + 1; ++Y) {
//
//					// Compute distance between object node and fluid lattice node.
//
//					const double dist_x = body_.at(i).cur_pos_.x_ - 0.5 - X;
//					const double dist_y = body_.at(i).cur_pos_.y_ + 0.5 - Y;
//
//					// Compute interpolation weights for x- and y-direction based on the distance.
//
//					const double weight_x = 1 - abs(dist_x);
//					const double weight_y = 1 - abs(dist_y);
//
//					// Compute lattice force.
//
//					fx(Y, (X + domain_x_) % domain_x_) += body_.at(i).Fx_ * weight_x * weight_y;
//					fy(Y, (X + domain_x_) % domain_x_) += body_.at(i).Fy_ * weight_x * weight_y;
//				}
//			}
//
//		}
//
//
//	}
//
//	void SpreadVelocity(Fluid & fluid)
//	{
//		for (int i = 0; i < nodes_num; ++i) 
//		{
//
//			// Reset node velocity first since '+=' is used.
//
//			body_.at(i).vx_ = 0.0;
//			body_.at(i).vy_ = 0.0;
//
//			// Identify the lowest fluid lattice node in interpolation range (see spreading).
//
//			int x_int = (int)(body_.at(i).cur_pos_.x_ - 0.5 + domain_x_) - domain_x_;
//			int y_int = (int)(body_.at(i).cur_pos_.y_ + 0.5);
//
//			// Run over all neighboring fluid nodes.
//			// In the case of the two-point interpolation, it is 2x2 fluid nodes.
//
//			for (int X = x_int; X <= x_int + 1; ++X) {
//				for (int Y = y_int; Y <= y_int + 1; ++Y) {
//
//					// Compute distance between object node and fluid lattice node.
//
//					const double dist_x = body_.at(i).cur_pos_.x_ - 0.5 - X;
//					const double dist_y = body_.at(i).cur_pos_.y_ + 0.5 - Y;
//
//					// Compute interpolation weights for x- and y-direction based on the distance.
//
//					const double weight_x = 1 - abs(dist_x);
//					const double weight_y = 1 - abs(dist_y);
//
//					// Compute node velocities.
//
//					body_.at(i).vx_ += (fluid.vx_(Y, (X + domain_x_) % domain_x_) * weight_x * weight_y);
//					body_.at(i).vy_ += (fluid.vy_(Y, (X + domain_x_) % domain_x_) * weight_x * weight_y);
//				}
//			}
//		}
//	}
//
//
//
//	void UpdatePosition()
//	{
//		/// Reset center position
//
//		center_.x_ = 0.0;
//		center_.y_ = 0.0;
//
//		/// Update node and center positions
//
//		for (int i = 0; i < nodes_num; ++i) 
//		{
//			body_.at(i).cur_pos_.x_ += body_.at(i).vx_;
//			body_.at(i).cur_pos_.y_ += body_.at(i).vy_;
//			
//			center_.x_ += body_.at(i).cur_pos_.x_ / nodes_num;
//			center_.y_ += body_.at(i).cur_pos_.y_ / nodes_num;
//		}
//
//		/// Check for periodicity along the x-axis
//
//		/*if (center_.x_ < 0) 
//		{
//			center_.x_ += domain_x_;
//
//			for (int n = 0; n < nodes_num; ++n) 
//			{
//				body_.at(n).cur_pos_.x_ += domain_x_;
//			}
//		}
//		else if (center_.x_ >= domain_x_) 
//		{
//			center_.x_ -= domain_x_;
//
//			for (int n = 0; n < nodes_num; ++n)
//			{
//				body_.at(n).cur_pos_.x_ -= domain_x_;
//			}
//		}*/
//
//
//	}
//
//	//! Writes data about boundary of immersed body to *.txt file
//	void WriteBodyFormToTxt(const int time)
//	{
//		std::string file_name = "Data/ib_lbm_data/body_form_txt/body_form_t" + std::to_string(time) + ".txt";
//
//		std::ofstream output_file;
//		output_file.open(file_name);
//
//		if (output_file.is_open())
//		{
//			for (auto node : body_)
//			{
//				output_file << node.cur_pos_.x_ << " " << node.cur_pos_.y_ << std::endl;
//			}
//
//			// Add first point to data in file to display closed boundary
//			output_file << body_.begin()->cur_pos_.x_ << " " << body_.begin()->cur_pos_.y_;
//
//		}
//		else
//			std::cout << "Error! Could not open file " << file_name << " to write form of immersed body. \n";
//
//		output_file.close();
//	}
//
//	//! Writes data about boundary of immersed body to *.vtk file
//	void WriteBodyFormToVtk(std::string file_path, const int time) 
//	{
//		std::string file_name = file_path + "\\body_fluid_t" + std::to_string(time) + ".vtk";
//
//		std::ofstream output_file;
//		output_file.open(file_name);
//
//		if (output_file.is_open())
//		{
//			// Write VTK header
//			output_file << "# vtk DataFile Version 3.0\n";
//			output_file << "particle_state\n";
//			output_file << "ASCII\n";
//			output_file << "DATASET POLYDATA\n";
//
//			// Write node positions
//			output_file << "POINTS " << nodes_num << " float\n";
//
//			for (auto node : body_)
//				output_file << node.cur_pos_.x_ << " " << node.cur_pos_.y_ << " 0\n";
//
//			// Write lines between neighboring nodes
//			output_file << "LINES " << nodes_num << " " << 3 * nodes_num << "\n";
//			
//			for (int i = 0; i < nodes_num; ++i) 
//				output_file << "2 " << i << " " << (i + 1) % nodes_num << "\n";
//
//			// Write vertices
//			output_file << "VERTICES 1 " << nodes_num + 1 << "\n";
//			output_file << nodes_num << " ";
//
//			for (int i = 0; i < nodes_num; ++i) 
//				output_file << i << " ";
//		}
//		else
//			std::cout << "Error! Could not open file " << file_name << " to write form of immersed body. \n";
//
//		output_file.close();
//	}
//
//private:
//
//	//! Calculate strain forces for all nodes of immersed body
//	void CalculateStrainForces()
//	{
//		for (int i = 0; i < nodes_num; ++i)
//		{
//			const double distance = SQ(body_.at(i).cur_pos_.x_ - body_.at((i + 1) % nodes_num).cur_pos_.x_) + SQ(body_.at(i).cur_pos_.y_ - body_.at((i + 1) % nodes_num).cur_pos_.y_);
//			const double distance_ref = SQ(body_.at(i).ref_pos_.x_ - body_.at((i + 1) % nodes_num).ref_pos_.x_) + SQ(body_.at(i).ref_pos_.y_ - body_.at((i + 1) % nodes_num).ref_pos_.y_);
//
//			const double fx = stiffness_ * (distance - distance_ref) * ( body_.at(i).cur_pos_.x_ - body_.at((i + 1) % nodes_num).cur_pos_.x_);
//			const double fy = stiffness_ * (distance - distance_ref) * ( body_.at(i).cur_pos_.y_ - body_.at((i + 1) % nodes_num).cur_pos_.y_);
//
//			// Signs of forces are chosen to satisfy third Newton law
//			body_.at(i).Fx_ += -fx;
//			body_.at(i).Fy_ += -fy;
//
//			body_.at((i + 1) % nodes_num).Fx_ += fx;
//			body_.at((i + 1) % nodes_num).Fy_ += fy;
//		}
//	}
//
//	void CalculateBendingForces()
//	{
//		for (int i = 0; i < nodes_num; ++i)
//		{
//			int prevId = (i - 1 + nodes_num) % nodes_num;
//			int nextId = (i + 1) % nodes_num;
//
//			const double x_l = body_.at(prevId).cur_pos_.x_;
//			const double y_l = body_.at(prevId).cur_pos_.y_;
//			const double x_m = body_.at(i).cur_pos_.x_;
//			const double y_m = body_.at(i).cur_pos_.y_;
//			const double x_r = body_.at(nextId).cur_pos_.x_;
//			const double y_r = body_.at(nextId).cur_pos_.y_;
//
//			const double x_l_ref = body_.at(prevId).ref_pos_.x_;
//			const double y_l_ref = body_.at(prevId).ref_pos_.y_;
//			const double x_m_ref = body_.at(i).ref_pos_.x_;
//			const double y_m_ref = body_.at(i).ref_pos_.y_;
//			const double x_r_ref = body_.at(nextId).ref_pos_.x_;
//			const double y_r_ref = body_.at(nextId).ref_pos_.y_;
//
//
//			// x-координата вектора, соединяющая l и r
//			const double tang_x_ref = x_r_ref - x_l_ref;
//			// y-координата вектора, соединяющая l и r
//			const double tang_y_ref = y_r_ref - y_l_ref;
//			double normal_x_ref;
//			double normal_y_ref;
//
//			// Тут просто задем нормаль так чтобы скалярное произведение вектора нормали на 
//			// вектор разности l и r были перпендикулярны : а разность модулей  для того чтобы 
//			// навпраление было всегда от вне
//			if (abs(tang_x_ref) < abs(tang_y_ref)) {
//				normal_x_ref = 1;
//				normal_y_ref = -tang_x_ref / tang_y_ref;
//			}
//			else {
//				normal_y_ref = 1;
//				normal_x_ref = -tang_y_ref / tang_x_ref;
//			}
//
//			// То же самое для обычныых
//			const double tang_x = x_r - x_l;
//			const double tang_y = y_r - y_l;
//			double normal_x;
//			double normal_y;
//
//			if (abs(tang_x) < abs(tang_y)) {
//				normal_x = 1;
//				normal_y = -tang_x / tang_y;
//			}
//			else {
//				normal_y = 1;
//				normal_x = -tang_y / tang_x;
//			}
//
//			// Просто нормализация вектора
//			const double normal_length_ref = sqrt(SQ(normal_x_ref) + SQ(normal_y_ref));
//			normal_x_ref /= normal_length_ref;
//			normal_y_ref /= normal_length_ref;
//
//			if (normal_x_ref * tang_y_ref - normal_y_ref * tang_x_ref > 0) {
//				normal_x_ref *= -1;
//				normal_y_ref *= -1;
//			}
//
//			const double normal_length = sqrt(SQ(normal_x) + SQ(normal_y));
//			normal_x /= normal_length;
//			normal_y /= normal_length;
//
//			if (normal_x * tang_y - normal_y * tang_x > 0) {
//				normal_x *= -1;
//				normal_y *= -1;
//			}
//			// Angle calculations
//
//			// Скалярное произведение векторов деленное на длину векторов
//			double angle_ref_cos = (x_l_ref - x_m_ref) * (x_m_ref - x_r_ref) + (y_l_ref - y_m_ref) * (y_m_ref - y_r_ref);
//			angle_ref_cos /= (sqrt(SQ(x_l_ref - x_m_ref) + SQ(y_l_ref - y_m_ref)) * sqrt(SQ(x_m_ref - x_r_ref) + SQ(y_m_ref - y_r_ref)));
//
//			// Addition if because of surrounding cos = 1.000000002
//			if (angle_ref_cos > 1.0)
//				angle_ref_cos = 1.0;
//			else if (angle_ref_cos < -1.0)
//				angle_ref_cos = -1.0;
//
//
//			double angle_ref = acos(angle_ref_cos);
//
//			const double convex_x_ref = (x_l_ref + x_r_ref) / 2 - x_m_ref;
//			const double convex_y_ref = (y_l_ref + y_r_ref) / 2 - y_m_ref;
//
//			if (convex_x_ref * normal_x_ref + convex_y_ref * normal_y_ref > 0) {
//				angle_ref *= -1;
//			}
//
//			double angle_cos = (x_l - x_m) * (x_m - x_r) + (y_l - y_m) * (y_m - y_r);
//			angle_cos /= (sqrt(SQ(x_l - x_m) + SQ(y_l - y_m)) * sqrt(SQ(x_m - x_r) + SQ(y_m - y_r)));
//
//			// Addition if because of surrounding cos = 1.000000002
//			if (angle_cos > 1.0)
//				angle_cos = 1.0;
//			else if (angle_cos < -1.0)
//				angle_cos = -1.0;
//
//
//			double angle = acos(angle_cos);
//
//			const double convex_x = (x_l + x_r) / 2 - x_m;
//			const double convex_y = (y_l + y_r) / 2 - y_m;
//
//			if (convex_x * normal_x + convex_y * normal_y > 0) {
//				angle *= -1;
//			}
//			// force calculations 
//			const double force_mag = bending_ * (angle - angle_ref);
//			const double length_l = abs(tang_x * (x_m - x_l) + tang_y * (y_m - y_l));
//			const double length_r = abs(tang_x * (x_m - x_r) + tang_y * (y_m - y_r));
//
//			body_.at(prevId).Fx_ += normal_x * force_mag * length_l / (length_l + length_r);
//			body_.at(prevId).Fy_ += normal_y * force_mag * length_l / (length_l + length_r);
//
//			body_.at(i).Fx_ += -normal_x * force_mag;
//			body_.at(i).Fy_ += -normal_y * force_mag;
//
//			body_.at(nextId).Fx_ += normal_x * force_mag * length_r / (length_l + length_r);
//			body_.at(nextId).Fy_ += normal_y * force_mag * length_r / (length_l + length_r);
//		}
//	}
//
//
//private:
//
//	int domain_x_;
//	int domain_y_;
//	//! Number of Lagragian nodes
//	int nodes_num;
//	//! Stiffness modulus 
//	const double stiffness_ = 0.1;
//	//! Bending modulus
//	const double bending_ = 0.001;
//
//	//! Radius of body
//	double radius_;
//	//! Center position of body
//	Point center_;
//
//	//! Nodes of Lagragian greed
//	std::vector<IBNode> body_;
//
//};



//class IBSolver
//{
//public:
//
//	IBSolver(double tau, Fluid& fluid, Medium & medium, ImmersedBody& body) : tau_(tau)
//	{ 
//		CreateDataFolder("Data\\ib_lbm_data");
//		CreateDataFolder("Data\\ib_lbm_data\\body_form_txt");
//		CreateDataFolder("Data\\ib_lbm_data\\body_form_vtk");
//		CreateDataFolder("Data\\ib_lbm_data\\fluid_txt");
//		CreateDataFolder("Data\\ib_lbm_data\\fluid_vtk");
//
//
//		std::cout << " --- Input parameters :\n";
//		std::cout << "nu = " << (tau - 0.5) / 3.0 << std::endl;
//
//		fluid_ = std::unique_ptr<Fluid>(new Fluid(fluid));
//		medium_ = std::unique_ptr<Medium>(new Medium(medium));
//		body_ = std::unique_ptr<ImmersedBody>(new ImmersedBody(body));
//
//		int rows = fluid_->size().first;
//		int colls = fluid_->size().second;
//		
//		fx_ = std::make_unique<Matrix2D<double>>(rows, colls);
//		fy_ = std::make_unique<Matrix2D<double>>(rows, colls);
//
//		force_member_.resize(kQ, 0.0);
//	}
//
//	~IBSolver() {}
//
//	//! Equilibrium probability distribution function calculation in SRT implementation
//	void feqCalculate()
//	{
//		// Проверить надо ли, или без нее все нормально
//		fluid_->feq_.fillWithoutBoundaries(0.0);
//
//		for (int q = 0; q < kQ; ++q) {
//			Matrix2D<double> v(fluid_->size().first, fluid_->size().second);
//			v = fluid_->vx_ * kEx[q] + fluid_->vy_ * kEy[q];
//
//			fluid_->feq_[q] = kW[q] * fluid_->rho_.ScalarMultiplication(
//				(1.0 + 3.0 * v + 4.5 * v.ScalarMultiplication(v) - 1.5 *
//				(fluid_->vx_.ScalarMultiplication(fluid_->vx_) + fluid_->vy_.ScalarMultiplication(fluid_->vy_)))
//			);
//		}
//	}
//
//	//! Streaming of particles to neighbour nodes in SRT implementation
//	void streaming()
//	{
//		for (int q = 0; q < kQ; ++q)
//		{
//			Matrix2D<double> temp = fluid_->f_[q];
//			fluid_->f_[q].FillWith(0.0);
//
//			for (unsigned y = 0; y < fluid_->size().first; ++y)
//				for (unsigned x = 0; x < fluid_->size().second; ++x)
//					if (medium_->is_fluid(y, x))
//						fluid_->f_[q](y - kEy[q], x + kEx[q]) = temp(y, x);
//		}
//
//		fluid_->f_.fillBoundaries(0.0);
//	}
//
//	void CalculateForces()
//	{
//		double gravity = 0.0;
//
//		for(int y = 0; y < fluid_->size().first; ++y)
//			for (int x = 0; x < fluid_->size().second; ++x)
//			{
//				force_member_.at(0) = (1.0 - 0.5 / tau_) * kW[0] * (3.0 * ((-fluid_->vx_(y, x) * ((*fx_)(y,x) + gravity) + -fluid_->vy_(y, x) * (*fy_)(y, x))));
//				force_member_[1] = (1 - 0.5 / tau_) * kW[1] * (3.0 * ((1 - fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity) + (-fluid_->vy_(y,x)) * (*fy_)(y,x)) + 9.0 * (fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity));
//				force_member_[2] = (1 - 0.5 / tau_) * kW[2] * (3.0* ((-1 - fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity) + (-fluid_->vy_(y,x)) * (*fy_)(y,x)) + 9.0 * (fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity));
//				force_member_[3] = (1 - 0.5 / tau_) * kW[3] * (3.0 * ((-fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity) + (1 - fluid_->vy_(y,x)) * (*fy_)(y,x)) + 9.0 * (fluid_->vy_(y,x)) * (*fy_)(y,x));
//				force_member_[4] = (1 - 0.5 / tau_) * kW[4] * (3.0 * ((-fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity) + (-1 - fluid_->vy_(y,x)) * (*fy_)(y,x)) + 9.0 * (fluid_->vy_(y,x)) * (*fy_)(y,x));
//				force_member_[5] = (1 - 0.5 / tau_) * kW[5] * (3.0 * ((1 - fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity) + (1 - fluid_->vy_(y,x)) * (*fy_)(y,x)) + 9.0 * (fluid_->vx_(y,x) + fluid_->vy_(y,x)) * ((*fx_)(y,x) + gravity + (*fy_)(y,x)));
//				force_member_[6] = (1 - 0.5 / tau_) * kW[6] * (3.0 * ((-1 - fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity) + (-1 - fluid_->vy_(y,x)) * (*fy_)(y,x)) + 9.0 * (fluid_->vx_(y,x) + fluid_->vy_(y,x)) * ((*fx_)(y,x) + gravity + (*fy_)(y,x)));
//				force_member_[7] = (1 - 0.5 / tau_) * kW[7] * (3.0 * ((1 - fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity) + (-1 - fluid_->vy_(y,x)) * (*fy_)(y,x)) + 9.0 * (fluid_->vx_(y,x) - fluid_->vy_(y,x)) * ((*fx_)(y,x) + gravity - (*fy_)(y,x)));
//				force_member_[8] = (1 - 0.5 / tau_) * kW[8] * (3.0 * ((-1 - fluid_->vx_(y,x)) * ((*fx_)(y,x) + gravity) + (1 - fluid_->vy_(y,x)) * (*fy_)(y,x)) + 9 * (fluid_->vx_(y,x) - fluid_->vy_(y,x)) * ((*fx_)(y,x) + gravity - (*fy_)(y,x)));
//			}
//	}
//
//	//! Collision of particles in nodes in SRT implementation
//	void collision()
//	{
//		CalculateForces();
//
//		for (int q = 0; q < kQ; ++q)
//			fluid_->f_[q] += (fluid_->feq_[q] - fluid_->f_[q]) / tau_ + force_member_.at(q);
//	}
//
//
//	//! Recalculation procedure (recalculate density, velocity) in SRT implementation
//	void recalculate()
//	{
//		fluid_->rho_ = fluid_->f_.calculateDensity();
//		fluid_->vx_ = fluid_->f_.calculateVelocity(kEx, fluid_->rho_, *fx_);
//
//		fluid_->vy_ = fluid_->f_.calculateVelocity(kEy, fluid_->rho_, *fy_);
//	}
//
//	void Solve(int iter_numb)
//	{
//		feqCalculate();
//
//		for (int q = 0; q < kQ; ++q)
//			fluid_->f_[q] = fluid_->feq_[q];
//
//		BCs BC(fluid_->f_);
//
//
//		for (int iter = 0; iter < iter_numb; ++iter)
//		{
//			body_->CalculateForces();
//			body_->SpreadForces(*fx_, *fy_);
//
//			collision();
//			BC.PrepareValuesForAllBC(BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::VON_NEUMAN, BCType::VON_NEUMAN);
//
//			streaming();
//
//			//BC.PrepareAdditionalBCs(*medium_);
//
//			BC.BounceBackBC(Boundary::TOP);
//			BC.BounceBackBC(Boundary::BOTTOM);
//			BC.VonNeumannBC(Boundary::LEFT, *fluid_, 0.01, 0.0);
//			BC.VonNeumannBC(Boundary::RIGHT, *fluid_, 0.01, 0.0);
//			//BC.BounceBackBC(Boundary::RIGHT);
//
//			//BC.AdditionalBounceBackBCs();
//
//			BC.RecordValuesForAllBC(BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::VON_NEUMAN, BCType::BOUNCE_BACK);
//
//			//BC.RecordAdditionalBCs();
//
//
//			recalculate();
//
//			feqCalculate();
//
//			body_->SpreadVelocity(*fluid_);
//			body_->UpdatePosition();
//
//			std::cout << iter << " Total rho = " << fluid_->rho_.GetSum() << std::endl;
//
//			if (iter % 25 == 0)
//			{
//				fluid_->write_fluid_vtk(iter);
//				body_->WriteBodyFormToTxt(iter);
//				body_->WriteBodyFormToVtk("Data\\ib_lbm_data\\body_form_vtk", iter);
//				fluid_->vx_.WriteFieldToTxt("Data\\ib_lbm_data\\fluid_txt", "vx", iter);
//				//fluid_->vx_.WriteToFile("vx", iter);
//			}
//
//		}
//
//		//! Streaming collision bounce back
//
//	}
//
//
//private:
//
//	//! Creates folder for output data if not existed yet
//	void CreateDataFolder(std::string folder_name) const
//	{
//		// Get path to current directory
//		char buffer[MAX_PATH];
//		GetModuleFileName(NULL, buffer, MAX_PATH);
//
//		std::string::size_type pos = std::string(buffer).find_last_of("\\/");
//		std::string path = std::string(buffer).substr(0, pos);
//		path = path.substr(0, path.size() - 6) + "\\" + folder_name;
//
//		char *cstr = new char[path.length() + 1];
//		strcpy(cstr, path.c_str());
//
//		// Create folder if not exist yet
//		if (GetFileAttributes(cstr) == INVALID_FILE_ATTRIBUTES) 
//			CreateDirectory(cstr, NULL);
//	}
//
//
//
//private:
//	//! Relaxation time
//	const double tau_;
//
//	std::unique_ptr<Fluid> fluid_;
//	std::unique_ptr<Medium> medium_;
//	std::unique_ptr<ImmersedBody> body_;
//
//	std::unique_ptr<Matrix2D<double>> fx_;
//	std::unique_ptr<Matrix2D<double>> fy_;
//
//	std::vector<double> force_member_;
//
//};

#pragma endregion

#pragma region 3d


// SRT approach implementation in 3D case.
//
// Relaxation parameter tau must be bigger then 0.5 to achive good results.
// It is better to choose it near 1.0;

class SRT3DSolver : iSolver
{
public:

	SRT3DSolver(double const tau, Medium3D & medium, Fluid3D & fluid);
	virtual ~SRT3DSolver() {}

	// Overriding iSolver methods
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



