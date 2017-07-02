#include<memory>
#include<windows.h>

#include <fstream> // file streams
#include <sstream> // string streams

#include"solver.h"
#include"..\modeling_area\fluid.h"
#include"..\modeling_area\medium.h"
#include"im_body\immersed_body.h"
#include"bc\bc.h"

#pragma region 2d


class IBSolver : iSolver
{
public:
	IBSolver(double tau, Fluid& fluid, Medium & medium, std::unique_ptr<ImmersedBody> body);
	IBSolver(double tau, Fluid& fluid, Medium & medium, std::vector<ImmersedBody*> bodies);

	//IBSolver(double tau, Fluid& fluid, Medium & medium, ImmersedBody& body);
	~IBSolver() {}

	void feqCalculate() override;
	void Streaming() override;

	void Collision() override;
	void Recalculate() override;
	void Solve(int iter_numb) override;

protected:

	//! Performs calculation of external force terms from immersed boundary on fluid
	void CalculateForces();
	//! Creates folder for output data if not existed yet
	void CreateDataFolder(std::string folder_name) const;

protected:
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

	// Add many immersed bodies
	std::vector<ImmersedBody*> im_bodies_;

};

class IBMRTSolver : IBSolver
{
public:

	IBMRTSolver(double tau, Fluid& fluid, Medium & medium, std::vector<ImmersedBody*> bodies) : IBSolver(tau, fluid, medium, bodies) 
	{
		CreateDataFolder("Data");
		CreateDataFolder("Data\\ib_lbm_data");
		CreateDataFolder("Data\\ib_lbm_data\\mrt");
		CreateDataFolder("Data\\ib_lbm_data\\mrt\\body_form_txt");
		CreateDataFolder("Data\\ib_lbm_data\\mrt\\body_form_vtk");
		CreateDataFolder("Data\\ib_lbm_data\\mrt\\fluid_txt");
		CreateDataFolder("Data\\ib_lbm_data\\mrt\\fluid_vtk");

		// Fill transformation matrix S
		S_ = { 1.0, 1.2, 1.0, 1.0, 1.2, 1.0, 1.2, 1.0 / tau_, 1.0 / tau_ };

		// Allocate memory for S * M^{-1} matrix
		MinvS_.Resize(kQ, kQ);

		// S - diagonal matrix
		Matrix2D<double> S(kQ, kQ);
		for (int i = 0; i < kQ; ++i)
			S(i, i) = S_.at(i);

		// Perform multiplication: S * M^{-1}
		for (int i = 0; i < kQ; ++i)
			for (int j = 0; j < kQ; ++j)
			{
				MinvS_(i, j) = 0.0;
				for (int k = 0; k < kQ; ++k)
					MinvS_(i, j) += Minv_[i][k] * S(k, j);
			}
	}

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
		dm[2] -= fluid_->rho_.ScalarMultiplication((-3.0 * vSq + 1));
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

	
	void Solve(int iter_numb) override
	{
		feqCalculate();

		for (int q = 0; q < kQ; ++q)
			fluid_->f_[q] = fluid_->feq_[q];

		BCs BC(fluid_->f_);

		for (int iter = 0; iter < iter_numb; ++iter)
		{
			// Clean fx, fy fields when many bodies are in modeling area [not in spread force function]
			fx_->FillWith(0.0);
			fy_->FillWith(0.0);

			for (auto& i : im_bodies_)
			{
				i->CalculateForces();
				//i->SpreadForces(*fx_, *fy_);
			}

			// Deal with RBC-Wall intearaction
			//Interaction(im_bodies_.at(2), im_bodies_.at(0));
			//Interaction(im_bodies_.at(2), im_bodies_.at(1));

			for (auto& i : im_bodies_)
			{
				//i->CalculateForces();
				i->SpreadForces(*fx_, *fy_);
			}

			Collision();
			BC.PrepareValuesForAllBC(BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::DIRICHLET, BCType::DIRICHLET);
			if (medium_->IsImmersedBodies())
				BC.PrepareAdditionalBCs(*medium_);

			Streaming();

			BC.BounceBackBC(Boundary::TOP);
			BC.BounceBackBC(Boundary::BOTTOM);

			BC.DirichletBC(Boundary::LEFT, *fluid_, 1.001);
			BC.DirichletBC(Boundary::RIGHT, *fluid_, 1.0);

			if (medium_->IsImmersedBodies())
				BC.AdditionalBounceBackBCs();

			BC.RecordValuesForAllBC(BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::DIRICHLET, BCType::DIRICHLET);
			if (medium_->IsImmersedBodies())
				BC.RecordAdditionalBCs();

			Recalculate();

			for (int y = 0; y < fluid_->GetRowsNumber(); ++y)
			{
				for (int x = 0; x < fluid_->GetColumnsNumber(); ++x)
					if (medium_->Get(y, x) == NodeType::OBSTACLE)
					{
						fluid_->vx_(y, x) = 0.0;
						fluid_->vy_(y, x) = 0.0;
						fluid_->rho_(y, x) = 0.0;
					}
			}

			feqCalculate();

			for (auto& i : im_bodies_)
			{
				i->SpreadVelocity(*fluid_);
				i->UpdatePosition();
			}

			std::cout << iter << " Total rho = " << fluid_->rho_.GetSum() << std::endl;

			if (iter % 50 == 0)
			{
				fluid_->WriteFluidToVTK("Data\\ib_lbm_data\\mrt\\fluid_vtk", iter);

				for (int i = 0; i < im_bodies_.size(); ++i)
				{
					im_bodies_.at(i)->WriteBodyFormToTxt(iter, i);
					im_bodies_.at(i)->WriteBodyFormToVtk("Data\\ib_lbm_data\\mrt\\body_form_vtk", i, iter);
				}

				fluid_->vx_.WriteFieldToTxt("Data\\ib_lbm_data\\mrt\\fluid_txt", "vx", iter);
				fluid_->vy_.WriteFieldToTxt("Data\\ib_lbm_data\\mrt\\fluid_txt", "vy", iter);
				fluid_->rho_.WriteFieldToTxt("Data\\ib_lbm_data\\mrt\\fluid_txt", "rho", iter);
			}

		}
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
