#include"srt.h"

SRTsolver::SRTsolver(double const tau, Medium & medium, Fluid & fluid) : 
	tau_(tau), 
	medium_(&medium), 
	fluid_(&fluid)
{
	assert(medium_->size().first == fluid_->size().first &&
		medium_->size().second == fluid_->size().second);
}

void SRTsolver::feq_calculate()
{
	// Проверить надо ли, или без нее все нормально
	fluid_->feq_.fill(0.0);

#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		Matrix<double> v(fluid_->size().first, fluid_->size().second);
		v = fluid_->vx_ * kEx[q] + fluid_->vy_ * kEy[q];

		fluid_->feq_[q] = kW[q] * fluid_->rho_.scalar_multiplication(
			(1.0 + 3.0 * v + 4.5 * v.scalar_multiplication(v) - 1.5 * 
			(fluid_->vx_.scalar_multiplication(fluid_->vx_) + fluid_->vy_.scalar_multiplication(fluid_->vy_)))
			);
	}
}

void SRTsolver::streaming()
{
	for (int q = 0; q < kQ; ++q) {
		Matrix<double> temp = fluid_->f_[q];
		fluid_->f_[q].fill_with(0.0);

		for (unsigned y = 0; y < fluid_->size().first; ++y)
			for (unsigned x = 0; x < fluid_->size().second; ++x)
				if (medium_->is_fluid(y, x))
			 		fluid_->f_[q](y - kEy[q], x + kEx[q]) = temp(y, x);
	}

	// Очищаем значения попавшие на границу, так как они уже сохранены в BCs
	fluid_->f_.fill_boundaries(0.0);
}

void SRTsolver::collision()
{
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] += (fluid_->feq_[q] - fluid_->f_[q]) / tau_;
}

void SRTsolver::solve(int iteration_number)
{
	Wall_info top{ Boundary::TOP, BCType::BOUNCE_BACK };
	Wall_info bottom{ Boundary::BOTTOM, BCType::BOUNCE_BACK };
	Wall_info left{ Boundary::LEFT, BCType::VON_NEUMAN };
	Wall_info right{ Boundary::RIGHT, BCType::BOUNCE_BACK };

	// check_bc()

	//std::cout << "initial rho = " << fluid_->rho_.get_sum() << std::endl;
	std::cout << fluid_->vx_;
	feq_calculate();

	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] = fluid_->feq_[q];

	BCs BC(fluid_->size().first, fluid_->size().second, fluid_->f_);

	for (int iter = 0; iter < iteration_number; ++iter) {
		collision();

		BC.prepair_bc_values(BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::VON_NEUMAN, BCType::BOUNCE_BACK);

		streaming();
		BC.bounce_back_bc(Boundary::TOP);
		BC.bounce_back_bc(Boundary::BOTTOM);
		BC.bounce_back_bc(Boundary::RIGHT);
		BC.von_neuman_bc(Boundary::LEFT, *fluid_, 0.01);


		BC.recording_bc_values(BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::VON_NEUMAN, BCType::BOUNCE_BACK);

		recalculate();

		std::vector<double> vel(fluid_->rows_ - 2, 0.01);
		fluid_->vx_.set_coll(1, vel);
		feq_calculate();

		std::cout << iter << " Total rho = " << fluid_->rho_.get_sum() << std::endl;
	}

	/*for (int i = 0; i < fluid_->rows_; ++i)
		std::cout << fluid_->vx_(i, 30) << "\n";*/
}

void SRTsolver::recalculate()
{
	fluid_->rho_ = fluid_->f_.get_density();
	fluid_->vx_ = fluid_->f_.get_velocity(kEx, fluid_->rho_);
	fluid_->vy_ = fluid_->f_.get_velocity(kEy, fluid_->rho_);
}


