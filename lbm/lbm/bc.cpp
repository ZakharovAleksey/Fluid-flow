#include"bc.h"

BCs::BCs(unsigned rows, unsigned colls, DistributionFunction<double> & dfunc): 
	length_(colls), height_(rows - 2), 
	f_ptr_(std::make_unique<DistributionFunction<double>>(dfunc)) {

	assert(length_ > 2 && height_ > 2);
#pragma omp parallel for
	for (int q = 0; q < kQ; ++q) {
		up_boundary_.at(q).resize(1, length_);
		bottom_boundary_.at(q).resize(1, length_);

		right_boundary_.at(q).resize(1, height_);
		left_boundary_.at(q).resize(1, height_);
	}

	std::cout << *f_ptr_;
}

BCs::~BCs() {}

bool BCs::get_values(Boundary BC)
{
	if (BC == Boundary::UP) {
		std::cout << "UP BOUNDARY VALUES:\n";
		bottom_boundary_ = f_ptr_->get_values_on_upper_boundary();
		return true;
	}
	else if (BC == Boundary::BOTTOM) {
		std::cout << "BOTTOM BOUNDARY VALUES:\n";
		bottom_boundary_ = f_ptr_->get_values_on_bottom_boundary();
		return true;
	}
	else if (BC == Boundary::LEFT) {
		std::cout << "LEFT BOUNDARY VALUES:\n";
		left_boundary_ = f_ptr_->get_values_on_left_boundary();
		return true;
	}
	else if (BC == Boundary::RIGHT) {
		std::cout << "RIGHT BOUNDARY VALUES:\n";
		right_boundary_ = f_ptr_->get_values_on_right_boundary();
		return true;
	}
	else
		return false;
}

bool BCs::get_boundady_values()
{
	return (get_values(Boundary::UP) && get_values(Boundary::BOTTOM) && 
		get_values(Boundary::LEFT) && get_values(Boundary::RIGHT));
}


void BCs::set_periodic_bc(Boundary first, Boundary second)
{
	if (first == Boundary::LEFT && second == Boundary::RIGHT) {
		right_boundary_.at(1).swap(left_boundary_.at(1));
		right_boundary_.at(5).swap(left_boundary_.at(5));
		right_boundary_.at(8).swap(left_boundary_.at(8));

		left_boundary_.at(3).swap(right_boundary_.at(3));
		left_boundary_.at(6).swap(right_boundary_.at(6));
		left_boundary_.at(7).swap(right_boundary_.at(7));
		
		std::cout << "PERIODIC BC LEFT RIGHT\n LEFT : \n";
		for (auto i : left_boundary_)
			std::cout << i;
		std::cout << "RIGHT\n";
		for (auto i : right_boundary_)
			std::cout << i;
	}
	else if (first == Boundary::UP && second == Boundary::BOTTOM) {
		bottom_boundary_.at(2).swap(bottom_boundary_.at(2));
		bottom_boundary_.at(5).swap(bottom_boundary_.at(5));
		bottom_boundary_.at(6).swap(bottom_boundary_.at(6));

		bottom_boundary_.at(4).swap(bottom_boundary_.at(4));
		bottom_boundary_.at(7).swap(bottom_boundary_.at(7));
		bottom_boundary_.at(8).swap(bottom_boundary_.at(8));
	
	}
	else {
		throw;
	}
}

void BCs::set_bounce_back_bc(Boundary first)
{
	if (first == Boundary::UP) {
		bottom_boundary_.at(4) = bottom_boundary_.at(2);
		bottom_boundary_.at(8) = bottom_boundary_.at(5);
		bottom_boundary_.at(7) = bottom_boundary_.at(6);

		// Пока что будем заполнять нулями для наглядности массивы, значения которых нам пока что не нужны

		bottom_boundary_.at(2).resize(1, 0);
		bottom_boundary_.at(5).resize(1, 0);
		bottom_boundary_.at(6).resize(1, 0);
		bottom_boundary_.at(1).resize(1, 0);
		bottom_boundary_.at(3).resize(1, 0);
		bottom_boundary_.at(0).resize(1, 0);

		/*std::cout << "BC BOUNCE BACK : UP BOUNDARY \n";
		for (int q = 0; q != kQ; ++q)
			std::cout << "f[" << q << "] = " << bottom_boundary_.at(q);*/
	}
	else if (first == Boundary::BOTTOM) {
		bottom_boundary_.at(2) = bottom_boundary_.at(4);
		bottom_boundary_.at(5) = bottom_boundary_.at(8);
		bottom_boundary_.at(6) = bottom_boundary_.at(7);

		// Пока что будем заполнять нулями для наглядности массивы, значения которых нам пока что не нужны

		bottom_boundary_.at(4).resize(1, 0);
		bottom_boundary_.at(7).resize(1, 0);
		bottom_boundary_.at(8).resize(1, 0);
		bottom_boundary_.at(1).resize(1, 0);
		bottom_boundary_.at(3).resize(1, 0);
		bottom_boundary_.at(0).resize(1, 0);

		/*std::cout << "BC BOUNCE BACK : BOTTOM BOUNDARY \n";
		for (int q = 0; q != kQ; ++q)
			std::cout << "f[" << q << "] = " << bottom_boundary_.at(q);*/
	}
	else if (first == Boundary::LEFT) {
		left_boundary_.at(1) = left_boundary_.at(3);
		left_boundary_.at(5) = left_boundary_.at(6);
		left_boundary_.at(8) = left_boundary_.at(7);

		// Пока что будем заполнять нулями для наглядности массивы, значения которых нам пока что не нужны
		left_boundary_.at(0).resize(1, 0);
		left_boundary_.at(2).resize(1, 0);
		left_boundary_.at(3).resize(1, 0);
		left_boundary_.at(4).resize(1, 0);
		left_boundary_.at(6).resize(1, 0);
		left_boundary_.at(7).resize(1, 0);

		/*std::cout << "BC BOUNCE BACK : LEFT BOUNDARY \n";
		for (int q = 0; q != kQ; ++q)
			std::cout << "f[" << q << "] = " << left_boundary_.at(q);*/
	}
	else if (first == Boundary::RIGHT) {

		right_boundary_.at(3) = right_boundary_.at(1);
		right_boundary_.at(6) = right_boundary_.at(5);
		right_boundary_.at(7) = right_boundary_.at(8);

		// Пока что будем заполнять нулями для наглядности массивы, значения которых нам пока что не нужны

		right_boundary_.at(0).resize(0, 0);
		std::cout << sizeof(right_boundary_.at(0)) << "---- \n";
		std::cout << right_boundary_.at(0).size().first << " sec = " << right_boundary_.at(0).size().second << std::endl;
		right_boundary_.at(1).resize(0, 0);
		right_boundary_.at(2).resize(0, 0);
		right_boundary_.at(4).resize(0, 0);
		right_boundary_.at(5).resize(0, 0);
		right_boundary_.at(8).resize(0, 0);

		std::cout << "BC BOUNCE BACK : RIGHT BOUNDARY \n";
		for (int q = 0; q != kQ; ++q)
			std::cout << "f[" << q << "] = " << right_boundary_.at(q);

	}
	else
		throw;

}

void BCs::set_von_neumann_bc(Boundary first, Fluid & fluid, double const v)
{
	if (first == Boundary::LEFT) {
		// Матрица в которую будут помещены знаения плотности
		Matrix<double> rho_temp(1, fluid.size().first - 2);

		rho_temp = (left_boundary_[0] + left_boundary_[2] + left_boundary_[4] +
			2.0 * (left_boundary_[3] + left_boundary_[6] + left_boundary_[7])) / (1.0 - v);

		left_boundary_[1] = left_boundary_[3] + (rho_temp * v * 2.0 / 3.0);

		Matrix<double> value = (left_boundary_[4] - left_boundary_[2]) / 2.0;

		left_boundary_[5] = rho_temp * v / 6.0 + left_boundary_[7] + value;
		left_boundary_[8] = rho_temp * v / 6.0 + left_boundary_[6] - value;

		fluid.rho_.set_coll(1, rho_temp);

		// Матрица которая хранит в себе все значения равные v чтобы заполнить скорости начальными знаениями
		Matrix<double> v_temp(1, fluid.size().first - 2); v_temp.fill_with(v);
		fluid.vx_.set_coll(1, v_temp);
	}
}

std::ostream & operator<<(std::ostream & os, BCs const & BC)
{
	os << "UP BOUNDARY ------ \n";
	for (int q = 0; q < kQ; ++q)
		os << BC.bottom_boundary_.at(q);

	os << "BOTTOM BOUNDARY ------ \n";
	for (int q = 0; q < kQ; ++q)
		os << BC.bottom_boundary_.at(q);

	os << "RIGHT BOUNDARY ------ \n";
	for (int q = 0; q < kQ; ++q)
		os << BC.right_boundary_.at(q);

	os << "LEFT BOUNDARY ------ \n";
	for (int q = 0; q < kQ; ++q)
		os << BC.left_boundary_.at(q);

	return os;
}
