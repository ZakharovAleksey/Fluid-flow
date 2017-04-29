#include"srt.h"

#pragma region 2d


SRTsolver::SRTsolver(double const tau, Medium & medium, Fluid & fluid) :
	tau_(tau),
	medium_(&medium),
	fluid_(&fluid)
{
	assert(medium_->size().first == fluid_->size().first &&
		medium_->size().second == fluid_->size().second);
}

void SRTsolver::feqCalculate()
{
	// Проверить надо ли, или без нее все нормально
	fluid_->feq_.fillWithoutBoundaries(0.0);

	for (int q = 0; q < kQ; ++q) {
		Matrix2D<double> v(fluid_->size().first, fluid_->size().second);
		v = fluid_->vx_ * kEx[q] + fluid_->vy_ * kEy[q];

		fluid_->feq_[q] = kW[q] * fluid_->rho_.ScalarMultiplication(
			(1.0 + 3.0 * v + 4.5 * v.ScalarMultiplication(v) - 1.5 *
			(fluid_->vx_.ScalarMultiplication(fluid_->vx_) + fluid_->vy_.ScalarMultiplication(fluid_->vy_)))
		);
	}
}

void SRTsolver::streaming()
{
	for (int q = 0; q < kQ; ++q) {
		Matrix2D<double> temp = fluid_->f_[q];
		fluid_->f_[q].FillWith(0.0);

		for (unsigned y = 0; y < fluid_->size().first; ++y)
			for (unsigned x = 0; x < fluid_->size().second; ++x)
				if (medium_->is_fluid(y, x))
					fluid_->f_[q](y - kEy[q], x + kEx[q]) = temp(y, x);
	}

	// Очищаем значения попавшие на границу, так как они уже сохранены в BCs
	fluid_->f_.fillBoundaries(0.0);
}

void SRTsolver::collision()
{
	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] += (fluid_->feq_[q] - fluid_->f_[q]) / tau_;
}

void SRTsolver::solve(int iteration_number)
{
	feqCalculate();
	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] = fluid_->feq_[q];

	BCs BC(fluid_->size().first, fluid_->size().second, fluid_->f_);

	for (int iter = 0; iter < iteration_number; ++iter) {
		collision();
		BC.prepareValuesForBC(BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::VON_NEUMAN, BCType::BOUNCE_BACK);

		streaming();

		BC.bounceBackBC(Boundary::TOP);
		BC.bounceBackBC(Boundary::BOTTOM);
		BC.bounceBackBC(Boundary::RIGHT);
		std::vector<double> vx;
		BC.vonNeumannBC(Boundary::LEFT, *fluid_, 0.01, vx);


		BC.recordValuesForBC(BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::VON_NEUMAN, BCType::BOUNCE_BACK);

		recalculate();
		fluid_->vx_.SetColumn(1, vx);

		feqCalculate();

		std::cout << iter << " Total rho = " << fluid_->rho_.GetSum() << std::endl;
	}

	for (int i = 0; i < fluid_->rows_; ++i)
		std::cout << fluid_->vx_(i, 5) << "\n";

	fluid_->vx_.WriteColumnToFile("vx", 5, 100);
	fluid_->vx_.WriteToFile("vx", iteration_number);
	fluid_->vx_.WriteRowToFile("vx", 5, 100);
}

void SRTsolver::recalculate()
{
	fluid_->rho_ = fluid_->f_.calculateDensity();
	fluid_->vx_ = fluid_->f_.calculateVelocity(kEx, fluid_->rho_);
	fluid_->vy_ = fluid_->f_.calculateVelocity(kEy, fluid_->rho_);
}




#pragma endregion

#pragma region 3d


SRT3DSolver::SRT3DSolver(double const tau, Medium3D & medium, Fluid3D & fluid) : tau_(tau), medium_(& medium), fluid_(& fluid)
{
	assert(medium_->GetDepthNumber() == fluid_->GetDepthNumber());
	assert(medium_->GetRowsNumber() == fluid_->GetRowsNumber());
	assert(medium_->GetColumnsNumber() == fluid_->GetColumnsNumber());
}

void SRT3DSolver::feqCalculate()
{
	int depth = medium_->GetDepthNumber();
	int rows = medium_->GetRowsNumber();
	int colls = medium_->GetColumnsNumber();
	
	// Obtain weights for feq calculation in accordance with Dmitry Biculov article
	std::vector<double> w;
	FillWeightsFor3D(w);

	for (int q = 0; q < kQ3d; ++q) 
	{
		// Article : Dmitry Biculov (e_{i}, v) form eq. (3) 
		Matrix3D<double> v(depth, rows, colls);
		v = *fluid_->vx_ * (double) ex[q] + *fluid_->vy_ * (double) ey[q] + *fluid_->vz_ * (double) ez[q];
		// Article : Dmitry Biculov (e_{i}, v)^{2} form eq. (3) 
		Matrix3D<double> v_sq(depth, rows, colls);
		v_sq = v.ScalarMultiplication(v);
		// Article : Dmitry Biculov (v, v)^{2} form eq. (3) 
		// !!! Copy becauce constant values : deal with it !!!
		Matrix3D<double> v_x = *fluid_->vx_;
		Matrix3D<double> v_y = *fluid_->vy_;
		Matrix3D<double> v_z = *fluid_->vz_;

		Matrix3D<double> v_2(depth, rows, colls);
		v_2 = v_x.ScalarMultiplication(v_x) + v_y.ScalarMultiplication(v_y) + v_z.ScalarMultiplication(v_z);


		fluid_->feq_->operator[](q) = w[q] * fluid_->rho_->ScalarMultiplication(1.0 + 3.0 * v + 4.5 * v_sq - 1.5 * v_2);
	}

}

void SRT3DSolver::streaming()
{
	const int depth = medium_->GetDepthNumber();
	const int rows = medium_->GetRowsNumber();
	const int colls = medium_->GetColumnsNumber();

	SubStreamingMiddle(depth, rows, colls);
	SubStreamingTop(depth, rows, colls);
	SubStreamingBottom(depth, rows, colls);
}

void SRT3DSolver::collision()
{
	for (int q = 0; q < kQ3d; ++q)
		(*fluid_->f_)[q] += ((*fluid_->feq_)[q] - (*fluid_->f_)[q]) / tau_;
}

void SRT3DSolver::solve(int iter_numb)
{
	fluid_->PoiseuilleIC(0.01);

	feqCalculate();
	for (int q = 0; q < kQ; ++q)
		(*fluid_->f_)[q] = (*fluid_->feq_)[q];

	BCs3D bc(fluid_->GetRowsNumber(), fluid_->GetColumnsNumber(), *fluid_->f_);

	for (int iter = 0; iter < iter_numb; ++iter)
	{
		std::cout << iter << " : ";
		collision();
		bc.PrepareValuesForBC(BCType::PERIODIC, BCType::PERIODIC, BCType::PERIODIC, BCType::PERIODIC, BCType::PERIODIC, BCType::PERIODIC);
		streaming();

		bc.PeriodicBC(Boundary::TOP, Boundary::BOTTOM);
		bc.PeriodicBC(Boundary::LEFT, Boundary::RIGHT);
		bc.PeriodicBC(Boundary::NEAR, Boundary::FAAR);

		bc.RecordValuesForBC(BCType::PERIODIC, BCType::PERIODIC, BCType::PERIODIC, BCType::PERIODIC, BCType::PERIODIC, BCType::PERIODIC);

		recalculate();
		feqCalculate();
	}
	GetProfile(5);
}

void SRT3DSolver::GetProfile(const int chan_numb)
{
	std::vector<double> res = fluid_->vy_->GetNFLayer(chan_numb);
	int colls = fluid_->GetColumnsNumber();
	int rows = fluid_->GetRowsNumber();
	int depth = fluid_->GetDepthNumber();

	std::ofstream file;
	file.open("Data/ex.txt");

	file.precision(3);
	//int i = 0;

	for (int i = 0; i < res.size(); ++i)
	{
		
		if(i != 0 && (i % (colls - 2) == 0))
			file << std::endl;
		file << res.at(i) << " ";
	}

	file.close();
	
	
}

void SRT3DSolver::SubStreamingMiddle(const int depth, const int rows, const int colls)
{
	for (int z = 0; z < depth; ++z)
	{
		for (int q = 0; q < 9; ++q)
		{
			Matrix2D<double> temp = fluid_->GetDistributionFuncLayer(z, q);
			fluid_->SetDistributionFuncLayerValue(z, q, 0.0);

			for (unsigned y = 0; y < rows; ++y)
				for (unsigned x = 0; x < colls; ++x)
					if (medium_->IsFluid(z, y, x))
						fluid_->f_->operator[](q)(z + ez[q], y + ey[q], x + ex[q]) = temp(y, x);
		}
	}
}

void SRT3DSolver::SubStreamingTop(const int depth, const int rows, const int colls)
{
	for (int z = 0; z < depth - 1; ++z)
	{
		for (int q = 9; q < 14; ++q)
		{
			Matrix2D<double> temp = fluid_->GetDistributionFuncLayer(z, q);
			fluid_->SetDistributionFuncLayerValue(z, q, 0.0);

			for (unsigned y = 0; y < rows; ++y)
				for (unsigned x = 0; x < colls; ++x)
					if (medium_->IsFluid(z, y, x))
						fluid_->f_->operator[](q)(z + ez[q], y + ey[q], x + ex[q]) = temp(y, x);
		}
	}
}

void SRT3DSolver::SubStreamingBottom(const int depth, const int rows, const int colls)
{
	for (int z = depth - 1; z > 0; --z)
	{
		for (int q = 14; q < 19; ++q)
		{
			Matrix2D<double> temp = fluid_->GetDistributionFuncLayer(z, q);
			fluid_->SetDistributionFuncLayerValue(z, q, 0.0);

			for (unsigned y = 0; y < rows; ++y)
				for (unsigned x = 0; x < colls; ++x)
					if (medium_->IsFluid(z, y, x))
						fluid_->f_->operator[](q)(z + ez[q], y + ey[q], x + ex[q]) = temp(y, x);
		}
	}
}

void SRT3DSolver::recalculate()
{
	fluid_->RecalculateRho();
	fluid_->RecalculateV();

	std::cout << "Total rho: " << fluid_->TotalRho() << std::endl;
}

#pragma endregion
