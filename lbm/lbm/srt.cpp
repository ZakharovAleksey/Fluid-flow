#include"srt.h"

SRTsolver::SRTsolver(double const tau, Medium & medium, Fluid & fluid) : 
	tau_(tau), 
	medium_(std::make_unique<Medium>(medium)), 
	fluid_(std::make_unique<Fluid>(fluid))
{
	assert(medium_->size().first == fluid_->size().first &&
		medium_->size().second == fluid_->size().second);
}

void SRTsolver::feq_calculate()
{
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
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		Matrix<double> temp = fluid_->f_[q];
		fluid_->f_[q].fill_with(0.0);

		for (unsigned y = 0; y < fluid_->size().first; ++y)
			for (unsigned x = 0; x < fluid_->size().second; ++x)
				if (medium_->is_fluid(y, x))
			 		fluid_->f_[q](y + kEy[q], x + kEx[q]) = temp(y, x);
	}
}

void SRTsolver::collision()
{
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] += (fluid_->feq_[q] - fluid_->f_[q]) / tau_;
}

void SRTsolver::solve()
{
	feq_calculate();
	collision();

	BCs BC(fluid_->size().first, fluid_->size().second, fluid_->f_);
	//if (BC.get_boundady_values()) {

	//	streaming();

	//	// --- Здесть задаются ГУ для каждой из стенок

	//	//BC.set_bounce_back_bc(Boundary::UP);
	//	//BC.set_bounce_back_bc(Boundary::BOTTOM);
	//	//BC.set_bounce_back_bc(Boundary::LEFT);
	//	//BC.set_bounce_back_bc(Boundary::RIGHT);
	//	
	//	//BC.set_von_neumann_bc(Boundary::LEFT, *fluid_, 0.01);
	//	// ---


	//}

	

}


