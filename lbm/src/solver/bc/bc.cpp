#include <functional>


#include"bc.h"
#include"../../math/array_func_impl.h" // Для работы с ГУ Фон-Неймана


#pragma region 2d




BCs::BCs(unsigned rows, unsigned colls, DistributionFunction<double> & dfunc): 
	length_(colls), height_(rows - 2), 
	f_ptr_(&dfunc) 
{
	assert(length_ > 2 && height_ > 2);
}

BCs::~BCs() {}

bool BCs::prepareValuesOnCurrentBoundary(Boundary const BC, BCType const boundary_condition_type)
{
	if (BC == Boundary::TOP) {
		if (boundary_condition_type == BCType::PERIODIC || boundary_condition_type == BCType::BOUNCE_BACK) {
			top_boundary_.insert(std::make_pair(2, f_ptr_->getTopBoundaryValues(2)));
			top_boundary_.insert(std::make_pair(5, f_ptr_->getTopBoundaryValues(5)));
			top_boundary_.insert(std::make_pair(6, f_ptr_->getTopBoundaryValues(6)));
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			// Пока что не реализовано исполькование ГУ типа Фон-Неймана
			throw;
		}

		return true;
	}
	else if (BC == Boundary::BOTTOM) {
		if (boundary_condition_type == BCType::PERIODIC || boundary_condition_type == BCType::BOUNCE_BACK) {
			bottom_boundary_.insert(std::make_pair(4, f_ptr_->getBottomBoundaryValue(4)));
			bottom_boundary_.insert(std::make_pair(7, f_ptr_->getBottomBoundaryValue(7)));
			bottom_boundary_.insert(std::make_pair(8, f_ptr_->getBottomBoundaryValue(8)));
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			// Пока что не реализовано исполькование ГУ типа Фон-Неймана
			throw;
		}

		return true;
	}
	else if (BC == Boundary::LEFT) {
		if (boundary_condition_type == BCType::PERIODIC || boundary_condition_type == BCType::BOUNCE_BACK) {
			left_boundary_.insert(std::make_pair(3, f_ptr_->getLeftBoundaryValue(3)));
			left_boundary_.insert(std::make_pair(6, f_ptr_->getLeftBoundaryValue(6)));
			left_boundary_.insert(std::make_pair(7, f_ptr_->getLeftBoundaryValue(7)));
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			left_boundary_.insert(std::make_pair(0, f_ptr_->getLeftBoundaryValue(0)));
			left_boundary_.insert(std::make_pair(2, f_ptr_->getLeftBoundaryValue(2)));
			left_boundary_.insert(std::make_pair(3, f_ptr_->getLeftBoundaryValue(3)));
			left_boundary_.insert(std::make_pair(4, f_ptr_->getLeftBoundaryValue(4)));
			left_boundary_.insert(std::make_pair(6, f_ptr_->getLeftBoundaryValue(6)));
			left_boundary_.insert(std::make_pair(7, f_ptr_->getLeftBoundaryValue(7)));
		}

		return true;
	}
	else if (BC == Boundary::RIGHT) {
		if (boundary_condition_type == BCType::PERIODIC || boundary_condition_type == BCType::BOUNCE_BACK) {
			right_boundary_.insert(std::make_pair(1, f_ptr_->getRightBoundaryValue(1)));
			right_boundary_.insert(std::make_pair(5, f_ptr_->getRightBoundaryValue(5)));
			right_boundary_.insert(std::make_pair(8, f_ptr_->getRightBoundaryValue(8)));
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

void BCs::prepareValuesForBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc)
{
	if (prepareValuesOnCurrentBoundary(Boundary::TOP, top_bc) &&
		prepareValuesOnCurrentBoundary(Boundary::BOTTOM, bottm_bc) &&
		prepareValuesOnCurrentBoundary(Boundary::LEFT, left_bc) &&
		prepareValuesOnCurrentBoundary(Boundary::RIGHT, right_bc)) 
	{
		// Лог что все значения получилось взять
	}
	else
	{
		throw;
	}
}

void BCs::recordValuesOnCurrentBoundary(Boundary const BC, BCType const boundary_condition_type)
{
	if (BC == Boundary::TOP) {
		if (boundary_condition_type == BCType::PERIODIC) {
			// auto = std::map<int, std::vector<double> >::iterator 
			auto f_2 = top_boundary_.find(2);
			auto f_5 = top_boundary_.find(5);
			auto f_6 = top_boundary_.find(6);

			f_ptr_->setBottomBoundaryValue(f_2->first, f_2->second);
			f_ptr_->setBottomBoundaryValue(f_5->first, f_5->second);
			f_ptr_->setBottomBoundaryValue(f_6->first, f_6->second);
		}
		else if (boundary_condition_type == BCType::BOUNCE_BACK) {
			auto f_4 = top_boundary_.find(4);
			auto f_7 = top_boundary_.find(7);
			auto f_8 = top_boundary_.find(8);

			f_ptr_->setTopBoundaryValue(f_4->first, f_4->second);
			f_ptr_->setTopBoundaryValue(f_7->first, f_7->second);
			f_ptr_->setTopBoundaryValue(f_8->first, f_8->second);
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

			f_ptr_->setTopBoundaryValue(f_4->first, f_4->second);
			f_ptr_->setTopBoundaryValue(f_7->first, f_7->second);
			f_ptr_->setTopBoundaryValue(f_8->first, f_8->second);
		}
		else if (boundary_condition_type == BCType::BOUNCE_BACK) {
			auto f_2 = bottom_boundary_.find(2);
			auto f_5 = bottom_boundary_.find(5);
			auto f_6 = bottom_boundary_.find(6);

			f_ptr_->setBottomBoundaryValue(f_2->first, f_2->second);
			f_ptr_->setBottomBoundaryValue(f_5->first, f_5->second);
			f_ptr_->setBottomBoundaryValue(f_6->first, f_6->second);
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

			f_ptr_->setRightBoundaryValue(f_3->first, f_3->second);
			f_ptr_->setRightBoundaryValue(f_6->first, f_6->second);
			f_ptr_->setRightBoundaryValue(f_7->first, f_7->second);
		}
		else if (boundary_condition_type == BCType::BOUNCE_BACK || 
				 boundary_condition_type == BCType::VON_NEUMAN) 
		{
			auto f_1 = left_boundary_.find(1);
			auto f_5 = left_boundary_.find(5);
			auto f_8 = left_boundary_.find(8);

			f_ptr_->setLeftBoundaryValue(f_1->first, f_1->second);
			f_ptr_->setLeftBoundaryValue(f_5->first, f_5->second);
			f_ptr_->setLeftBoundaryValue(f_8->first, f_8->second);
		}
	}
	else if (BC == Boundary::RIGHT) {
		if (boundary_condition_type == BCType::PERIODIC) {
			auto f_1 = right_boundary_.find(1);
			auto f_5 = right_boundary_.find(5);
			auto f_8 = right_boundary_.find(8);

			f_ptr_->setLeftBoundaryValue(f_1->first, f_1->second);
			f_ptr_->setLeftBoundaryValue(f_5->first, f_5->second);
			f_ptr_->setLeftBoundaryValue(f_8->first, f_8->second);
		}
		else if (boundary_condition_type == BCType::BOUNCE_BACK) {
			auto f_3 = right_boundary_.find(3);
			auto f_6 = right_boundary_.find(6);
			auto f_7 = right_boundary_.find(7);

			f_ptr_->setRightBoundaryValue(f_3->first, f_3->second);
			f_ptr_->setRightBoundaryValue(f_6->first, f_6->second);
			f_ptr_->setRightBoundaryValue(f_7->first, f_7->second);
		}
		else if (boundary_condition_type == BCType::VON_NEUMAN) {
			// Пока еще не реализованно
			throw;
		}
	}
}

void BCs::recordValuesForBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc)
{
	recordValuesOnCurrentBoundary(Boundary::TOP, top_bc);
	recordValuesOnCurrentBoundary(Boundary::BOTTOM, bottm_bc);
	recordValuesOnCurrentBoundary(Boundary::LEFT, left_bc);
	recordValuesOnCurrentBoundary(Boundary::RIGHT, right_bc);
}

void BCs::periodicBC(Boundary const first, Boundary const second)
{
	if (first == Boundary::LEFT && second == Boundary::RIGHT)
		left_boundary_.swap(right_boundary_);
	else if (first == Boundary::TOP && second == Boundary::BOTTOM)
		top_boundary_.swap(bottom_boundary_);
	else
		throw;
}

void BCs::bounceBackBC(Boundary const first)
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

void BCs::vonNeumannBC(Boundary const first, Fluid & fluid, double const vx, 
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


#pragma endregion


BCs3D::BCs3D(int rows, int colls, DistributionFunction3D<double>& dfunc) :
	height_(rows), length_(colls - 2), f_ptr_(&dfunc)
{
	
}

bool BCs3D::PrepareValuesInCurBoundary(Boundary const BC, BCType const bc_type)
{
	// Pointer to function, which gets appropriate values, bepending on wall type:
	// for TOP : ptr to GetTopBoundaryValues, for BOTTOM ptr to GetBottomBoundaryValue
	std::vector<double>(DistributionFunction3D<double>::*ptrToFunc)(int) const = nullptr;


	switch (BC)
	{
	case Boundary::TOP:

		ptrToFunc = &DistributionFunction3D<double>::GetTopBoundaryValues;
		return WriteBoundaryValues(bc_type, top_boundary_, top_ids_, ptrToFunc);
		break;

	case Boundary::BOTTOM:

		ptrToFunc = &DistributionFunction3D<double>::GetBottomBoundaryValue;
		return WriteBoundaryValues(bc_type, bottom_boundary_, bottom_ids_, ptrToFunc);
		break;

	case Boundary::RIGHT:

		ptrToFunc = &DistributionFunction3D<double>::GetRightBoundaryValue;
		return WriteBoundaryValues(bc_type, right_boundary_, right_ids_, ptrToFunc);

		break;

	case Boundary::LEFT:

		ptrToFunc = &DistributionFunction3D<double>::GetLeftBoundaryValue;
		return WriteBoundaryValues(bc_type, left_boundary_, left_ids_, ptrToFunc);
		break;

	case Boundary::NEAR:

		ptrToFunc = &DistributionFunction3D<double>::GetNearBoundaryValue;
		return WriteBoundaryValues(bc_type, near_boundary_, near_ids_, ptrToFunc);
		break;

	case Boundary::FAAR:

		ptrToFunc = &DistributionFunction3D<double>::GetFarBoundaryValue;
		return WriteBoundaryValues(bc_type, far_boundary_, far_ids_, ptrToFunc);
		break;

	default:
		std::cout << "Wrong Boundary type is used while prepair values for BC.\n";
		return false;
		break;
	}
}

bool BCs3D::WriteBoundaryValues(BCType const bc_type, std::map<int, std::vector<double> > & bc_boundary, const std::vector<int> & bc_ids, 
	/* Pointer to function to distinguish boundaries: TOP, BOTTOM, e.t.c */std::vector<double> (DistributionFunction3D<double>::*ptrToFunc)(int) const)
{
	if (bc_type == BCType::PERIODIC || bc_type == BCType::BOUNCE_BACK)
	{
		for (auto bc_id : bc_ids)
			// Get appropriate values from probability distribution function (using APPROPRIATE function via ptr. to function) 
			// and write them in appropriate class fields.
			bc_boundary.insert(std::make_pair(bc_id, (*f_ptr_.*ptrToFunc)(bc_id)));

		return true;

	}
	else if (bc_type == BCType::VON_NEUMAN)
	{
		std::cout << "Von Neumann BC are not implemented yet!\n";

		return false;
		throw;
	}
	else
	{
		std::cout << "Wrond BC Type is used!\n";
	}
	return false;
}

void BCs3D::PrepareValuesForBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc, BCType const near_bc, BCType far_bc)
{
	if (PrepareValuesInCurBoundary(Boundary::TOP, top_bc) &&
		PrepareValuesInCurBoundary(Boundary::BOTTOM, bottm_bc) &&
		PrepareValuesInCurBoundary(Boundary::LEFT, left_bc) &&
		PrepareValuesInCurBoundary(Boundary::RIGHT, right_bc) && 
		PrepareValuesInCurBoundary(Boundary::NEAR, near_bc) &&
		PrepareValuesInCurBoundary(Boundary::FAAR, far_bc) )
	{
		// Лог что все значения получилось взять
	}
	else
	{
		std::cout << "Could not prepare values for all boundaries.\n";
		throw;
	}
}

void BCs3D::recordValuesOnCurrentBoundary(Boundary const BC, BCType const bc_type)
{
	// Берем значения с одной стороны -а засовывем в другую! в периодических берем с топ - суем в ботом
	if (BC == Boundary::TOP) 
	{
		if (bc_type == BCType::PERIODIC) 
		{
			// auto = std::map<int, std::vector<double> >::iterator 

			for (auto topId : top_ids_)
			{
				auto f_id = bottom_boundary_.find(topId);
				f_ptr_->SetBottomBoundaryValue(f_id->first, f_id->second);
			}
			/*auto f_9 = bottom_boundary_.find(9);
			auto f_10 = bottom_boundary_.find(10);
			auto f_11 = bottom_boundary_.find(11);
			auto f_12 = bottom_boundary_.find(12);
			auto f_13 = bottom_boundary_.find(13);

			f_ptr_->SetBottomBoundaryValue(f_9->first, f_9->second);
			f_ptr_->SetBottomBoundaryValue(f_10->first, f_10->second);
			f_ptr_->SetBottomBoundaryValue(f_11->first, f_11->second);
			f_ptr_->SetBottomBoundaryValue(f_12->first, f_12->second);
			f_ptr_->SetBottomBoundaryValue(f_13->first, f_13->second);*/
		}
		else if (bc_type == BCType::BOUNCE_BACK) 
		{
			// !!! ОСТАНОВИЛСЯ НА СОЗДАНИИ УНИВЕРСАЛЬНОЙ ФУНКЦИИ ДЛЯ BB УСЛОВИЙ ЧЕРЕЗ УКАЗАТЕЛЬ НА ФУНКЦИЮ!!!
			for (auto botId : bottom_ids_)
			{
				auto f_id = top_boundary_.find(botId);
				f_ptr_->SetTopBoundaryValue(f_id->first, f_id->second);
			}
		}
		else if (bc_type == BCType::VON_NEUMAN) 
		{
			std::cout << "Von NEUMAN on TOP is not yet implemented\n";
			throw;
		}
	}
	
	else if (BC == Boundary::BOTTOM) 
	{
		if (bc_type == BCType::PERIODIC) 
		{
			auto f_14 = top_boundary_.find(14);
			auto f_15 = top_boundary_.find(15);
			auto f_16 = top_boundary_.find(16);
			auto f_17 = top_boundary_.find(17);
			auto f_18 = top_boundary_.find(18);

			f_ptr_->SetTopBoundaryValue(f_14->first, f_14->second);
			f_ptr_->SetTopBoundaryValue(f_15->first, f_15->second);
			f_ptr_->SetTopBoundaryValue(f_16->first, f_16->second);
			f_ptr_->SetTopBoundaryValue(f_17->first, f_17->second);
			f_ptr_->SetTopBoundaryValue(f_18->first, f_18->second);
		}
		else if (bc_type == BCType::BOUNCE_BACK) 
		{
			for (auto topId : top_ids_)
			{
				auto f_id = bottom_boundary_.find(topId);
				f_ptr_->SetBottomBoundaryValue(f_id->first, f_id->second);
			}
		}
		else if (bc_type == BCType::VON_NEUMAN) {
			// Пока еще не реализованно
			throw;
		}
	}
	
	else if (BC == Boundary::RIGHT) 
	{
		if (bc_type == BCType::PERIODIC) {
			auto f_1 = left_boundary_.find(1);
			auto f_5 = left_boundary_.find(5);
			auto f_8 = left_boundary_.find(8);
			auto f_10 = left_boundary_.find(10);
			auto f_15 = left_boundary_.find(15);

			f_ptr_->SetLeftBoundaryValue(f_1->first, f_1->second);
			f_ptr_->SetLeftBoundaryValue(f_5->first, f_5->second);
			f_ptr_->SetLeftBoundaryValue(f_8->first, f_8->second);
			f_ptr_->SetLeftBoundaryValue(f_10->first, f_10->second);
			f_ptr_->SetLeftBoundaryValue(f_15->first, f_15->second);
		}
		else if (bc_type == BCType::BOUNCE_BACK)
		{
			for (auto leftId : left_ids_)
			{
				auto f_id = right_boundary_.find(leftId);
				f_ptr_->SetRightBoundaryValue(f_id->first, f_id->second);
			}
		}
		else if (bc_type == BCType::VON_NEUMAN)
		{
			throw;
		}
	}
	
	else if (BC == Boundary::LEFT) 
	{
		if (bc_type == BCType::PERIODIC) 
		{
			auto f_3 = right_boundary_.find(3);
			auto f_6 = right_boundary_.find(6);
			auto f_7 = right_boundary_.find(7);
			auto f_12 = right_boundary_.find(12);
			auto f_17 = right_boundary_.find(17);

			f_ptr_->SetRightBoundaryValue(f_3->first, f_3->second);
			f_ptr_->SetRightBoundaryValue(f_6->first, f_6->second);
			f_ptr_->SetRightBoundaryValue(f_7->first, f_7->second);
			f_ptr_->SetRightBoundaryValue(f_12->first, f_12->second);
			f_ptr_->SetRightBoundaryValue(f_17->first, f_17->second);
		}
		else if (bc_type == BCType::BOUNCE_BACK) 
		{
			for (auto rightId : right_ids_)
			{
				auto f_id = left_boundary_.find(rightId);
				f_ptr_->SetLeftBoundaryValue(f_id->first, f_id->second);
			}
		}
		else if (bc_type == BCType::VON_NEUMAN) 
		{
			// Пока еще не реализованно
			throw;
		}
	}
	
	else if (BC == Boundary::NEAR) {
		if (bc_type == BCType::PERIODIC) {
			auto f_4 = far_boundary_.find(4);
			auto f_7 = far_boundary_.find(7);
			auto f_8 = far_boundary_.find(8);
			auto f_13 = far_boundary_.find(13);
			auto f_18 = far_boundary_.find(18);

			f_ptr_->SetFarBoundaryValue(f_4->first, f_4->second);
			f_ptr_->SetFarBoundaryValue(f_7->first, f_7->second);
			f_ptr_->SetFarBoundaryValue(f_8->first, f_8->second);
			f_ptr_->SetFarBoundaryValue(f_13->first, f_13->second);
			f_ptr_->SetFarBoundaryValue(f_18->first, f_18->second);
		}
		else if (bc_type == BCType::BOUNCE_BACK) 
		{
			for (auto farId : far_ids_)
			{
				auto f_id = near_boundary_.find(farId);
				f_ptr_->SetNearBoundaryValue(f_id->first, f_id->second);
			}
		}
		else if (bc_type == BCType::VON_NEUMAN) 
		{
			// Пока еще не реализованно
			throw;
		}
	}
	
	else if (BC == Boundary::FAAR) {
		if (bc_type == BCType::PERIODIC) 
		{
			auto f_2 = near_boundary_.find(2);
			auto f_5 = near_boundary_.find(5);
			auto f_6 = near_boundary_.find(6);
			auto f_11 = near_boundary_.find(11);
			auto f_16 = near_boundary_.find(16);

			f_ptr_->SetNearBoundaryValue(f_2->first, f_2->second);
			f_ptr_->SetNearBoundaryValue(f_5->first, f_5->second);
			f_ptr_->SetNearBoundaryValue(f_6->first, f_6->second);
			f_ptr_->SetNearBoundaryValue(f_11->first, f_11->second);
			f_ptr_->SetNearBoundaryValue(f_16->first, f_16->second);
		}
		else if (bc_type == BCType::BOUNCE_BACK) 
		{
			for (auto nearId : near_ids_)
			{
				auto f_id = far_boundary_.find(nearId);
				f_ptr_->SetFarBoundaryValue(f_id->first, f_id->second);
			}
		}
		else if (bc_type == BCType::VON_NEUMAN) 
		{
			// Пока еще не реализованно
			throw;
		}
	}
}

void BCs3D::RecordValuesForBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc, BCType const near_bc, BCType far_bc)
{
	recordValuesOnCurrentBoundary(Boundary::TOP, top_bc);
	recordValuesOnCurrentBoundary(Boundary::BOTTOM, bottm_bc);
	recordValuesOnCurrentBoundary(Boundary::LEFT, left_bc);
	recordValuesOnCurrentBoundary(Boundary::RIGHT, right_bc);
	recordValuesOnCurrentBoundary(Boundary::NEAR, near_bc);
	recordValuesOnCurrentBoundary(Boundary::FAAR, far_bc);

	f_ptr_->ClearBoundaries();
}

void BCs3D::PeriodicBC(Boundary const first, Boundary const second)
{
	if (first == Boundary::LEFT && second == Boundary::RIGHT)
		left_boundary_.swap(right_boundary_);
	else if (first == Boundary::TOP && second == Boundary::BOTTOM)
		top_boundary_.swap(bottom_boundary_);
	else if (first == Boundary::NEAR && second == Boundary::FAAR)
		near_boundary_.swap(far_boundary_);
	else
	{
		std::cout << "Check parameters in PeriodicBC function.\n";
		throw;
	}
		
}

void BCs3D::BounceBackBC(Boundary const first)
{
	if (first == Boundary::TOP) 
	{
		swap_id(top_boundary_, 9, 14);
		swap_id(top_boundary_, 10, 17);
		swap_id(top_boundary_, 11, 18);
		swap_id(top_boundary_, 12, 15);
		swap_id(top_boundary_, 13, 16);
	}
	else if (first == Boundary::BOTTOM) 
	{
		swap_id(bottom_boundary_, 14, 9);
		swap_id(bottom_boundary_, 15, 12);
		swap_id(bottom_boundary_, 16, 13);
		swap_id(bottom_boundary_, 17, 10);
		swap_id(bottom_boundary_, 18, 11);
	}
	else if (first == Boundary::LEFT) 
	{
		swap_id(left_boundary_, 3, 1);
		swap_id(left_boundary_, 6, 8);
		swap_id(left_boundary_, 7, 5);
		swap_id(left_boundary_, 12, 15);
		swap_id(left_boundary_, 17, 10);
	}
	else if (first == Boundary::RIGHT) 
	{
		swap_id(right_boundary_, 1, 3);
		swap_id(right_boundary_, 5, 7);
		swap_id(right_boundary_, 8, 6);
		swap_id(right_boundary_, 10, 17);
		swap_id(right_boundary_, 15, 12);
	}
	else if (first == Boundary::NEAR)
	{
		swap_id(near_boundary_, 4, 2);
		swap_id(near_boundary_, 7, 5);
		swap_id(near_boundary_, 8, 6);
		swap_id(near_boundary_, 13, 16);
		swap_id(near_boundary_, 18, 11);
	}
	else if (first == Boundary::FAAR)
	{
		swap_id(far_boundary_, 2, 4);
		swap_id(far_boundary_, 5, 7);
		swap_id(far_boundary_, 6, 8);
		swap_id(far_boundary_, 11, 18);
		swap_id(far_boundary_, 16, 13);
	}
}

void BCs3D::swap_id(std::map<int, std::vector<double>>& map, int const from, int const to)
{
	std::vector<double> temp;
	auto iter = map.find(from);
	temp.swap(iter->second);
	map.erase(iter);
	map.insert(std::make_pair(to, temp));
}
