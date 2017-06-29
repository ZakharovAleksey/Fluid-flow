#include"mrt.h"

MRTSolver::MRTSolver(double const tau, Medium & medium, Fluid & fluid) : SRTsolver(tau, medium, fluid)
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
	for (int i = 0; i < kQ; ++i)
		for (int j = 0; j < kQ; ++j)
		{
			MinvS_(i, j) = 0.0;
			for (int k = 0; k < kQ; ++k)
				MinvS_(i, j) += Minv_[i][k] * S(k, j);
		}

}

void MRTSolver::Collision()
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

void MRTSolver::Solve(int iteration_number)
{
	feqCalculate();
	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] = fluid_->feq_[q];

	BCs BC(fluid_->f_);

	BloodFlowMicrophone mic(Point(40, 100));

	for (int i = 0; i < 6; ++i)
	{
		double x = 35 + 20 * sin(i * M_PI / ( 2 * 5 ));
		double y = 39 - 20 * cos(i * M_PI / (2 * 5));
		mic.AddMeasurePoint(Point(y, x));
	}

	for (int i = 1; i < 6; ++i)
	{
		double x = 50 + 25 * i;
		double y = 35;
		mic.AddMeasurePoint(Point(y, x));
		y = 30;
		mic.AddMeasurePoint(Point(y, x));
	}

	for (int iter = 0; iter < iteration_number; ++iter)
	{
		mic.TakeOffParameters(iter, fluid_->vx_, "vx");
		mic.TakeOffParameters(iter, fluid_->rho_, "rho");

		Collision();
		BC.PrepareValuesForAllBC(BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::DIRICHLET, BCType::DIRICHLET);
		if (medium_->IsImmersedBodies())
			BC.PrepareAdditionalBCs(*medium_);

		Streaming();

		//BC.PeriodicBC(Boundary::TOP, Boundary::BOTTOM);
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
			for(int x = 0; x < fluid_->GetColumnsNumber(); ++x)
				if (medium_->Get(y, x) == NodeType::OBSTACLE)
				{
					fluid_->vx_(y, x) = 0.0;
					fluid_->vy_(y, x) = 0.0;
					fluid_->rho_(y, x) = 0.0;
				}
		}

		feqCalculate();

		std::cout << iter << " Total rho = " << fluid_->rho_.GetSum() << std::endl;

		if (iter % 100 == 0)
		{
			//Matrix2D<double> v = CalculateModulus(fluid_->vx_, fluid_->vy_);
			//v.WriteFieldToTxt("Data\\mrt_lbm_data\\2d\\fluid_txt", "v", iter);
			fluid_->vx_.WriteFieldToTxt("Data\\mrt_lbm_data\\2d\\fluid_txt", "vx", iter);
			fluid_->vy_.WriteFieldToTxt("Data\\mrt_lbm_data\\2d\\fluid_txt", "vy", iter);
			fluid_->WriteFluidToVTK("Data\\mrt_lbm_data\\2d\\fluid_vtk", iter);
		}

	}

}