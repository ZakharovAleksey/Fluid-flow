#include <functional>

#include"bc.h"

#include"array_func_impl.h" // Для работы с ГУ Фон-Неймана


BCs::BCs(unsigned rows, unsigned colls, DistributionFunction<double> & dfunc): 
	length_(colls), height_(rows - 2), 
	f_ptr_(&dfunc) 
{
	assert(length_ > 2 && height_ > 2);
}

BCs::~BCs() {}

bool BCs::get_values(Boundary const BC, BCType const boundary_condition_type)
{
	if (BC == Boundary::TOP) {
		if (boundary_condition_type == BCType::PERIODIC || boundary_condition_type == BCType::BOUNCE_BACK) {
			top_boundary_.insert(std::make_pair(2, f_ptr_->get_top_boundary_val(2)));
			top_boundary_.insert(std::make_pair(5, f_ptr_->get_top_boundary_val(5)));
			top_boundary_.insert(std::make_pair(6, f_ptr_->get_top_boundary_val(6)));
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			// Пока что не реализовано исполькование ГУ типа Фон-Неймана
			throw;
		}

		return true;
	}
	else if (BC == Boundary::BOTTOM) {
		if (boundary_condition_type == BCType::PERIODIC || boundary_condition_type == BCType::BOUNCE_BACK) {
			bottom_boundary_.insert(std::make_pair(4, f_ptr_->get_bottom_boundary_val(4)));
			bottom_boundary_.insert(std::make_pair(7, f_ptr_->get_bottom_boundary_val(7)));
			bottom_boundary_.insert(std::make_pair(8, f_ptr_->get_bottom_boundary_val(8)));
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			// Пока что не реализовано исполькование ГУ типа Фон-Неймана
			throw;
		}

		return true;
	}
	else if (BC == Boundary::LEFT) {
		if (boundary_condition_type == BCType::PERIODIC || boundary_condition_type == BCType::BOUNCE_BACK) {
			left_boundary_.insert(std::make_pair(3, f_ptr_->get_left_boundary_val(3)));
			left_boundary_.insert(std::make_pair(6, f_ptr_->get_left_boundary_val(6)));
			left_boundary_.insert(std::make_pair(7, f_ptr_->get_left_boundary_val(7)));
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			left_boundary_.insert(std::make_pair(0, f_ptr_->get_left_boundary_val(0)));
			left_boundary_.insert(std::make_pair(2, f_ptr_->get_left_boundary_val(2)));
			left_boundary_.insert(std::make_pair(3, f_ptr_->get_left_boundary_val(3)));
			left_boundary_.insert(std::make_pair(4, f_ptr_->get_left_boundary_val(4)));
			left_boundary_.insert(std::make_pair(6, f_ptr_->get_left_boundary_val(6)));
			left_boundary_.insert(std::make_pair(7, f_ptr_->get_left_boundary_val(7)));
		}

		return true;
	}
	else if (BC == Boundary::RIGHT) {
		if (boundary_condition_type == BCType::PERIODIC || boundary_condition_type == BCType::BOUNCE_BACK) {
			right_boundary_.insert(std::make_pair(1, f_ptr_->get_right_boundary_val(1)));
			right_boundary_.insert(std::make_pair(5, f_ptr_->get_right_boundary_val(5)));
			right_boundary_.insert(std::make_pair(8, f_ptr_->get_right_boundary_val(8)));
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			// Пока что не реализовано исполькование ГУ типа Фон-Неймана
			throw;
		}

		return true;
	}
	else
		return false;
}

void BCs::prepair_bc_values(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc)
{
	if (get_values(Boundary::TOP, top_bc) &&
		get_values(Boundary::BOTTOM, bottm_bc) &&
		get_values(Boundary::LEFT, left_bc) &&
		get_values(Boundary::RIGHT, right_bc)) 
	{
		// Лог что все значения получилось взять
	}
	else
	{
		throw;
	}
}

void BCs::set_values(Boundary const BC, BCType const boundary_condition_type)
{
	if (BC == Boundary::TOP) {
		if (boundary_condition_type == BCType::PERIODIC) {
			// auto = std::map<int, std::vector<double> >::iterator 
			auto f_2 = top_boundary_.find(2);
			auto f_5 = top_boundary_.find(5);
			auto f_6 = top_boundary_.find(6);

			f_ptr_->set_bottom_boundary_value(f_2->first, f_2->second);
			f_ptr_->set_bottom_boundary_value(f_5->first, f_5->second);
			f_ptr_->set_bottom_boundary_value(f_6->first, f_6->second);
		}
		else if (boundary_condition_type == BCType::BOUNCE_BACK) {
			auto f_4 = top_boundary_.find(4);
			auto f_7 = top_boundary_.find(7);
			auto f_8 = top_boundary_.find(8);

			f_ptr_->set_top_boundary_value(f_4->first, f_4->second);
			f_ptr_->set_top_boundary_value(f_7->first, f_7->second);
			f_ptr_->set_top_boundary_value(f_8->first, f_8->second);
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			// Пока еще не реализованно
			throw;
		}
	}
	else if (BC == Boundary::BOTTOM) {
		if (boundary_condition_type == BCType::PERIODIC) {
			auto f_4 = bottom_boundary_.find(4);
			auto f_7 = bottom_boundary_.find(7);
			auto f_8 = bottom_boundary_.find(8);

			f_ptr_->set_top_boundary_value(f_4->first, f_4->second);
			f_ptr_->set_top_boundary_value(f_7->first, f_7->second);
			f_ptr_->set_top_boundary_value(f_8->first, f_8->second);
		}
		else if (boundary_condition_type == BCType::BOUNCE_BACK) {
			auto f_2 = bottom_boundary_.find(2);
			auto f_5 = bottom_boundary_.find(5);
			auto f_6 = bottom_boundary_.find(6);

			f_ptr_->set_bottom_boundary_value(f_2->first, f_2->second);
			f_ptr_->set_bottom_boundary_value(f_5->first, f_5->second);
			f_ptr_->set_bottom_boundary_value(f_6->first, f_6->second);
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			// Пока еще не реализованно
			throw;
		}
	}
	else if (BC == Boundary::LEFT) {
		if (boundary_condition_type == BCType::PERIODIC) {
			auto f_3 = left_boundary_.find(3);
			auto f_6 = left_boundary_.find(6);
			auto f_7 = left_boundary_.find(7);

			f_ptr_->set_right_boundary_value(f_3->first, f_3->second);
			f_ptr_->set_right_boundary_value(f_6->first, f_6->second);
			f_ptr_->set_right_boundary_value(f_7->first, f_7->second);
		}
		else if (boundary_condition_type == BCType::BOUNCE_BACK || 
				 boundary_condition_type == BCType::VON_NEUMAN) 
		{
			auto f_1 = left_boundary_.find(1);
			auto f_5 = left_boundary_.find(5);
			auto f_8 = left_boundary_.find(8);

			f_ptr_->set_left_boundary_value(f_1->first, f_1->second);
			f_ptr_->set_left_boundary_value(f_5->first, f_5->second);
			f_ptr_->set_left_boundary_value(f_8->first, f_8->second);
		}
	}
	else if (BC == Boundary::RIGHT) {
		if (boundary_condition_type == BCType::PERIODIC) {
			auto f_1 = right_boundary_.find(1);
			auto f_5 = right_boundary_.find(5);
			auto f_8 = right_boundary_.find(8);

			f_ptr_->set_left_boundary_value(f_1->first, f_1->second);
			f_ptr_->set_left_boundary_value(f_5->first, f_5->second);
			f_ptr_->set_left_boundary_value(f_8->first, f_8->second);
		}
		else if (boundary_condition_type == BCType::BOUNCE_BACK) {
			auto f_3 = right_boundary_.find(3);
			auto f_6 = right_boundary_.find(6);
			auto f_7 = right_boundary_.find(7);

			f_ptr_->set_right_boundary_value(f_3->first, f_3->second);
			f_ptr_->set_right_boundary_value(f_6->first, f_6->second);
			f_ptr_->set_right_boundary_value(f_7->first, f_7->second);
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			// Пока еще не реализованно
			throw;
		}
	}
}

void BCs::recording_bc_values(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc)
{
	set_values(Boundary::TOP, top_bc);
	set_values(Boundary::BOTTOM, bottm_bc);
	set_values(Boundary::LEFT, left_bc);
	set_values(Boundary::RIGHT, right_bc);
}

void BCs::periodic_bc(Boundary const first, Boundary const second)
{
	if (first == Boundary::LEFT && second == Boundary::RIGHT)
		left_boundary_.swap(right_boundary_);
	else if (first == Boundary::TOP && second == Boundary::BOTTOM)
		top_boundary_.swap(bottom_boundary_);
	else
		throw;
}

void BCs::bounce_back_bc(Boundary const first)
{
	if (first == Boundary::TOP) {
		swap_id(top_boundary_, 2, 4);
		swap_id(top_boundary_, 5, 8);
		swap_id(top_boundary_, 6, 7);
	}
	else if(first == Boundary::BOTTOM) {
		swap_id(bottom_boundary_, 4, 2);
		swap_id(bottom_boundary_, 8, 5);
		swap_id(bottom_boundary_, 7, 6);
	}
	else if (first == Boundary::LEFT) {
		swap_id(left_boundary_, 3, 1);
		swap_id(left_boundary_, 6, 8);
		swap_id(left_boundary_, 7, 5);
	}
	else if (first == Boundary::RIGHT) {
		swap_id(right_boundary_, 1, 3);
		swap_id(right_boundary_, 5, 7);
		swap_id(right_boundary_, 8, 6);
	}
}

void BCs::von_neuman_bc(Boundary const first, Fluid & fluid, double const vx, 
	std::vector<double> & velocity_x)
{
	// Подготовка векторов, куда запишутся скорость и плотность на границе
	if (vx != 0.0) {
		if (velocity_x.empty())
			velocity_x.resize(fluid.size().first - 2, vx);
		else {
			velocity_x.clear();
			velocity_x.resize(fluid.size().first - 2, vx);
		}
	}

	if (first == Boundary::LEFT) {
		std::vector<double> density(fluid.size().first - 2, 0.0);
		density = left_boundary_.at(0) + left_boundary_.at(2) + left_boundary_.at(4) +
			(left_boundary_.at(3) + left_boundary_.at(6) + left_boundary_.at(7)) * 2.0 / (1.0 - vx);

		left_boundary_.insert(std::make_pair(1, left_boundary_.at(3) + (density * vx * 2.0 / 3.0)));

		std::vector<double> value = (left_boundary_.at(4) - left_boundary_.at(2)) / 2.0;

		left_boundary_.insert(std::make_pair(5, density * vx / 6.0 + left_boundary_.at(7) + value));
		left_boundary_.insert(std::make_pair(8, density * vx / 6.0 + left_boundary_.at(6) - value));
		

		// Матрица которая хранит в себе все значения равные v чтобы заполнить скорости начальными знаениями

		// Теперь, когда вычесленны f[1], f[5], f[8], rho, v - удаляем ненужные поля
		left_boundary_.erase(0);
		left_boundary_.erase(2);
		left_boundary_.erase(3);
		left_boundary_.erase(4);
		left_boundary_.erase(6);
		left_boundary_.erase(7);
	}
}


void BCs::swap_id(std::map<int, std::vector<double>> & map, int const from, int const to)
{
	std::vector<double> temp;
	auto iter = map.find(from);
	temp.swap(iter->second);
	map.erase(iter);
	map.insert(std::make_pair(to, temp));

}

std::ostream & operator<<(std::ostream & os, BCs const & BC)
{
	os.precision(3);

	os << "TOP BOUNDARY ------ \n";
	for (auto i : BC.top_boundary_) {
		os << "f[" << i.first << "] = ";
		for (auto j : i.second)
			os << j << ' ';
		os << std::endl;
	}

	os << "BOTTOM BOUNDARY ------ \n";
	for (auto i : BC.bottom_boundary_) {
		os << "f[" << i.first << "] = ";
		for (auto j : i.second)
			os << j << ' ';
		os << std::endl;
	}

	os << "RIGHT BOUNDARY ------ \n";
	for (auto i : BC.right_boundary_) {
		os << "f[" << i.first << "] = ";
		for (auto j : i.second)
			os << j << ' ';
		os << std::endl;
	}

	os << "LEFT BOUNDARY ------ \n";
	for (auto i : BC.left_boundary_) {
		os << "f[" << i.first << "] = ";
		for (auto j : i.second)
			os << j << ' ';
		os << std::endl;
	}

	return os;
}
