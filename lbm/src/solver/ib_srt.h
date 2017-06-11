#include<memory>
#include<windows.h>

#include <fstream> // file streams
#include <sstream> // string streams

#include"solver.h"
#include"..\modeling_area\fluid.h"
#include"..\modeling_area\medium.h"
#include"im_body\immersed_body.h"
#include"bc\bc.h"

//#define 	M_PI   3.14159265358979323846
//#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code


class IBSolver : iSolver
{
public:
	IBSolver(double tau, Fluid& fluid, Medium & medium, std::unique_ptr<ImmersedBody> body);
	//IBSolver(double tau, Fluid& fluid, Medium & medium, ImmersedBody& body);
	~IBSolver() {}

	void feqCalculate() override;
	void Streaming() override;

	void Collision() override;
	void Recalculate() override;
	void Solve(int iter_numb) override;

private:

	//! Performs calculation of external force terms from immersed boundary on fluid
	void CalculateForces();

	//! Creates folder for output data if not existed yet
	void CreateDataFolder(std::string folder_name) const;

private:
	//! Relaxation time
	const double tau_;

	//! Pointer to fluid class of appropriate modeling area
	std::unique_ptr<Fluid> fluid_;
	//! Pointer to medium class of appropriate modeling area
	std::unique_ptr<Medium> medium_;
	//! Pointer to immersed body
	std::unique_ptr<ImmersedBody> body_;
	
	//! Pointer to x-component external forces for all modeling area
	std::unique_ptr<Matrix2D<double>> fx_;
	//! Pointer to y-component external forces for all modeling area
	std::unique_ptr<Matrix2D<double>> fy_;

	std::vector<double> force_member_;

};