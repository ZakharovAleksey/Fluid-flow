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

	for (int iter = 0; iter < iteration_number; ++iter)
	{
		Collision();
		BC.PrepareValuesForAllBC(BCType::VON_NEUMAN, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK);

		Streaming();

		BC.VonNeumannBC(Boundary::TOP, *fluid_, 0.01, 0.0);
		BC.BounceBackBC(Boundary::BOTTOM);
		BC.BounceBackBC(Boundary::LEFT);
		BC.BounceBackBC(Boundary::RIGHT);

		BC.RecordValuesForAllBC(BCType::VON_NEUMAN, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK);

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