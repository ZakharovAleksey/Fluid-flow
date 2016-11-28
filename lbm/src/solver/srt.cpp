#include"srt.h"

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
		Matrix<double> v(fluid_->size().first, fluid_->size().second);
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
		Matrix<double> temp = fluid_->f_[q];
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


