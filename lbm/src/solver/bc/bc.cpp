#include <functional>


#include"bc.h"
#include"../../math/array_func_impl.h" // Для работы с ГУ Фон-Неймана


#pragma region 2d


BCs::BCs(unsigned rows, unsigned colls, DistributionFunction<double> & dfunc): length_(colls), height_(rows - 2), f_ptr_(&dfunc) { }

BCs::~BCs() {}

bool BCs::PrepareValuesForSingleBC(Boundary const BC, BCType const bc_type)
{
	std::vector<double>(DistributionFunction<double>::*ptrToFunc)(int) const = nullptr;

	switch (BC)
	{
	case Boundary::TOP:

		ptrToFunc = &DistributionFunction<double>::getTopBoundaryValues;

		if (bc_type == BCType::PERIODIC || bc_type == BCType::BOUNCE_BACK)
			return WriteBoundaryValues(bc_type, top_boundary_, top_ids_, ptrToFunc);
		// !! IMPLEMENTED BUT NOT TESTED !!!
		else if (bc_type == BCType::VON_NEUMAN)
			return WriteVonNeumannBoundaryValues(bc_type, top_boundary_, top_ids_, mid_width_ids_, ptrToFunc);
		
		break;

	case Boundary::BOTTOM:

		ptrToFunc = &DistributionFunction<double>::getBottomBoundaryValue;

		if (bc_type == BCType::PERIODIC || bc_type == BCType::BOUNCE_BACK)
			return WriteBoundaryValues(bc_type, bottom_boundary_, bottom_ids_, ptrToFunc);
		// !! IMPLEMENTED BUT NOT TESTED !!!
		else if (bc_type == BCType::VON_NEUMAN)
			return WriteVonNeumannBoundaryValues(bc_type, bottom_boundary_, bottom_ids_, mid_width_ids_, ptrToFunc);

		break;

	case Boundary::LEFT:

		ptrToFunc = &DistributionFunction<double>::getLeftBoundaryValue;

		if (bc_type == BCType::PERIODIC || bc_type == BCType::BOUNCE_BACK)
			return WriteBoundaryValues(bc_type, left_boundary_, left_ids_, ptrToFunc);
		// !! IMPLEMENTED BUT NOT TESTED !!!
		else if (bc_type == BCType::VON_NEUMAN)
			return WriteVonNeumannBoundaryValues(bc_type, left_boundary_, left_ids_, mid_height_ids_, ptrToFunc);

		break;

	case Boundary::RIGHT:

		ptrToFunc = &DistributionFunction<double>::getRightBoundaryValue;

		if (bc_type == BCType::PERIODIC || bc_type == BCType::BOUNCE_BACK)
			return WriteBoundaryValues(bc_type, right_boundary_, right_ids_, ptrToFunc);
		// !! IMPLEMENTED BUT NOT TESTED !!!
		else if (bc_type == BCType::VON_NEUMAN)
			return WriteVonNeumannBoundaryValues(bc_type, right_boundary_, right_ids_, mid_height_ids_, ptrToFunc);

		break;

	default:
		std::cout << "Try to prepare values for wrong BC Type.\n";
		break;
	}
}

void BCs::PrepareValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc)
{
	if (PrepareValuesForSingleBC(Boundary::TOP, top_bc) &&
		PrepareValuesForSingleBC(Boundary::BOTTOM, bottm_bc) &&
		PrepareValuesForSingleBC(Boundary::LEFT, left_bc) &&
		PrepareValuesForSingleBC(Boundary::RIGHT, right_bc)) 
	{
		// Лог что все значения получилось взять
	}
	else
	{
		std::cout << "Error while prepairing BC values!\n";
		throw;
	}
}

void BCs::RecordValuesOnSingleBC(Boundary const BC, BCType const bc_type)
{
	void(DistributionFunction<double>::*ptrToFunc)(int const, std::vector<double> const &) = nullptr;

	switch (BC)
	{
	case Boundary::TOP:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction<double>::setBottomBoundaryValue;
			RecordBoundaryValues(bc_type, top_boundary_, top_ids_, ptrToFunc);
		}
		else if (bc_type == BCType::BOUNCE_BACK || bc_type == BCType::VON_NEUMAN) // !! VON NEUMANN IS NOT TESTED YET !!!
		{
			ptrToFunc = &DistributionFunction<double>::setTopBoundaryValue;
			RecordBoundaryValues(bc_type, top_boundary_, bottom_ids_, ptrToFunc);
		}

		break;
	case Boundary::BOTTOM:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction<double>::setTopBoundaryValue;
			RecordBoundaryValues(bc_type, bottom_boundary_, bottom_ids_, ptrToFunc);
		}
		else if (bc_type == BCType::BOUNCE_BACK || bc_type == BCType::VON_NEUMAN) // !! VON NEUMANN IS NOT TESTED YET !!!
		{
			ptrToFunc = &DistributionFunction<double>::setBottomBoundaryValue;
			RecordBoundaryValues(bc_type, bottom_boundary_, top_ids_, ptrToFunc);
		}

		break;
	case Boundary::LEFT:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction<double>::setRightBoundaryValue;
			RecordBoundaryValues(bc_type, left_boundary_, left_ids_, ptrToFunc);
		}
		else if (bc_type == BCType::BOUNCE_BACK || bc_type == BCType::VON_NEUMAN)
		{
			ptrToFunc = &DistributionFunction<double>::setLeftBoundaryValue;
			RecordBoundaryValues(bc_type, left_boundary_, right_ids_, ptrToFunc);
		}

		break;
	case Boundary::RIGHT:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction<double>::setLeftBoundaryValue;
			RecordBoundaryValues(bc_type, right_boundary_, right_ids_, ptrToFunc);
		}
		else if (bc_type == BCType::BOUNCE_BACK || bc_type == BCType::VON_NEUMAN) // !! VON NEUMANN IS NOT TESTED YET !!!
		{
			ptrToFunc = &DistributionFunction<double>::setRightBoundaryValue;
			RecordBoundaryValues(bc_type, right_boundary_, left_ids_, ptrToFunc);
		}

		break;
	default:
		std::cout << "Wrong BC type appears while record BC.\n";
		break;
	}
}

void BCs::RecordValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc)
{
	RecordValuesOnSingleBC(Boundary::TOP, top_bc);
	RecordValuesOnSingleBC(Boundary::BOTTOM, bottm_bc);
	RecordValuesOnSingleBC(Boundary::LEFT, left_bc);
	RecordValuesOnSingleBC(Boundary::RIGHT, right_bc);
}

bool BCs::WriteBoundaryValues(BCType const bc_type, std::map<int, std::vector<double>>& bc_boundary, const std::vector<int>& bc_ids, std::vector<double>(DistributionFunction<double>::* ptrToFunc)(int) const)
{
	if (bc_type == BCType::PERIODIC || bc_type == BCType::BOUNCE_BACK)
	{
		for (auto bc_id : bc_ids)
			// Get appropriate values from probability distribution function (using APPROPRIATE function via ptr. to function) 
			// and write them in appropriate class fields.
			bc_boundary.insert(std::make_pair(bc_id, (*f_ptr_.*ptrToFunc)(bc_id)));

		return true;

	}
	return false;
}

bool BCs::WriteVonNeumannBoundaryValues(BCType const bc_type, std::map<int, std::vector<double>>& bc_boundary, const std::vector<int>& bc_ids_1, const std::vector<int>& bc_ids_2, std::vector<double>(DistributionFunction<double>::* ptrToFunc)(int) const)
{
	if (bc_type == BCType::VON_NEUMAN)
	{
		// Idea: Get Values from first and second layers: for example for BOTTOM we need ALL distribution function from MIDDLE and TOP layers
		// Get appropriate values from probability distribution function (using APPROPRIATE function via ptr. to function) 
		// and write them in appropriate class fields.

		for (auto bc_id_1 : bc_ids_1)
			bc_boundary.insert(std::make_pair(bc_id_1, (*f_ptr_.*ptrToFunc)(bc_id_1)));
		for (auto bc_id_2 : bc_ids_2)
			bc_boundary.insert(std::make_pair(bc_id_2, (*f_ptr_.*ptrToFunc)(bc_id_2)));

		return true;
	}
	else
	{

	}
	return false;
}

bool BCs::RecordBoundaryValues(BCType const bc_type, std::map<int, std::vector<double>>& bc_boundary, const std::vector<int>& bc_ids, void(DistributionFunction<double>::* ptrToFunc)(int const, std::vector<double> const &))
{

	if (bc_type == BCType::PERIODIC || bc_type == BCType::BOUNCE_BACK) // !!! Так как одиноковые то может сюда добавить и bc_type == BCType::VON_NEUMAN !!!
	{
		for (auto bc_id : bc_ids)
		{
			auto f_id = bc_boundary.find(bc_id);
			(*f_ptr_.*ptrToFunc)(f_id->first, f_id->second);
		}
		return true;
	}
	else if (bc_type == BCType::VON_NEUMAN)
	{
		// !!! In testing mode
		for (auto bc_id : bc_ids)
		{
			auto f_id = bc_boundary.find(bc_id);
			(*f_ptr_.*ptrToFunc)(f_id->first, f_id->second);
		}
		return true;
	}
	else
	{
		std::cout << "Wrong BC Type is used!\n";
	}
	return false;
}


void BCs::PeriodicBC(Boundary const first, Boundary const second)
{
	if (first == Boundary::LEFT && second == Boundary::RIGHT)
		left_boundary_.swap(right_boundary_);
	else if (first == Boundary::TOP && second == Boundary::BOTTOM)
		top_boundary_.swap(bottom_boundary_);
	else
		throw;
}

void BCs::BounceBackBC(Boundary const first)
{
	if (first == Boundary::TOP) 
	{
		SwapId(top_boundary_, 2, 4);
		SwapId(top_boundary_, 5, 8);
		SwapId(top_boundary_, 6, 7);
	}
	else if(first == Boundary::BOTTOM) 
	{
		SwapId(bottom_boundary_, 4, 2);
		SwapId(bottom_boundary_, 8, 5);
		SwapId(bottom_boundary_, 7, 6);
	}
	else if (first == Boundary::LEFT) 
	{
		SwapId(left_boundary_, 3, 1);
		SwapId(left_boundary_, 6, 8);
		SwapId(left_boundary_, 7, 5);
	}
	else if (first == Boundary::RIGHT) 
	{
		SwapId(right_boundary_, 1, 3);
		SwapId(right_boundary_, 5, 7);
		SwapId(right_boundary_, 8, 6);
	}
}

void BCs::VonNeumannBC_OLD(Boundary const first, Fluid & fluid, double const vx, std::vector<double> & velocity_x)
{
	// Подготовка векторов, куда запишутся скорость и плотность на границе
	if (vx != 0.0) 
	{
		if (velocity_x.empty())
			velocity_x.resize(fluid.size().first - 2, vx);
		else {
			velocity_x.clear();
			velocity_x.resize(fluid.size().first - 2, vx);
		}
	}

	if (first == Boundary::LEFT) 
	{
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

void BCs::VonNeumannBC(Boundary const first, Fluid & fluid, double const vx, double const vy)
{
	switch (first)
	{
	case Boundary::TOP:
		VonNeumannBC1(first, fluid, mid_width_ids_, top_ids_, vx, vy);
		break;
	case Boundary::BOTTOM:
		VonNeumannBC1(first, fluid, mid_width_ids_, bottom_ids_, vx, vy);
		break;
	case Boundary::LEFT:
		VonNeumannBC1(first, fluid, mid_height_ids_, left_ids_, vx, vy);
		break;
	case Boundary::RIGHT:
		VonNeumannBC1(first, fluid, mid_height_ids_, right_ids_, vx, vy);
		break;
	default:
		std::cout << "Try to apply Von-Neumann BC to wrong boundary type.\n";
		break;
	}
}

void BCs::VonNeumannBC1(Boundary const first, Fluid & fluid, const std::vector<int> ids_1, const std::vector<int> ids_2, double const vx, double const vy)
{
	// Arrays with incoming velocities for further calcuations

	const int x_size = fluid.size().second;
	const int y_size = fluid.size().first - 2;

	std::vector<double> vel_x(x_size, vx);
	std::vector<double> vel_y(y_size, vy);

	if (first == Boundary::TOP)
	{
		std::vector<double> rho(x_size, 0.0);

		for (auto id : ids_1)
			rho = rho + top_boundary_.at(id);
		for(auto id : ids_2)
			rho = rho + top_boundary_.at(id) * 2.0;

		rho = rho / (1.0 + vy);

		top_boundary_.insert(std::make_pair(4, top_boundary_.at(2) - 2.0 / 3.0 * rho * vy));


		std::vector<double> temp = (top_boundary_.at(1) - top_boundary_.at(3)) / 2.0;

		top_boundary_.insert(std::make_pair(7, top_boundary_.at(5) + temp - (vy / 6.0 + vx / 2.0) * rho));
		top_boundary_.insert(std::make_pair(8, top_boundary_.at(6) - temp - (vy / 6.0 - vx / 2.0) * rho));

		// + set velocities to 1 row after recalculation
		for (auto id : ids_1)
			top_boundary_.erase(id);

		for (auto id : ids_2)
			top_boundary_.erase(id);
	}
	else if (first == Boundary::BOTTOM)
	{
		std::vector<double> rho(x_size, 0.0);

		for (auto id : ids_1)
			rho = rho + bottom_boundary_.at(id);
		for (auto id : ids_2)
			rho = rho + bottom_boundary_.at(id) * 2.0;

		rho = rho / (1.0 - vy);

		bottom_boundary_.insert(std::make_pair(2, bottom_boundary_.at(4) + 2.0 / 3.0 * rho * vy));


		std::vector<double> temp = (bottom_boundary_.at(1) - bottom_boundary_.at(3)) / 2.0;

		bottom_boundary_.insert(std::make_pair(5, bottom_boundary_.at(7) - temp + (vy / 6.0 + vx / 2.0) * rho));
		bottom_boundary_.insert(std::make_pair(6, bottom_boundary_.at(8) + temp + (vy / 6.0 - vx / 2.0) * rho));

		// + set velocities to fluid->size().second - 1 row after recalculation
		for (auto id : ids_1)
			bottom_boundary_.erase(id);

		for (auto id : ids_2)
			bottom_boundary_.erase(id);
	}
	else if (first == Boundary::LEFT)
	{
		std::vector<double> rho(y_size, 0.0);

		for (auto id : ids_1)
			rho = rho + left_boundary_.at(id);
		for (auto id : ids_2)
			rho = rho + left_boundary_.at(id) * 2.0;

		rho = rho / (1.0 - vx);

		left_boundary_.insert(std::make_pair(1, left_boundary_.at(3) + 2.0 / 3.0 * rho * vx));


		std::vector<double> temp = (left_boundary_.at(2) - left_boundary_.at(4)) / 2.0;

		left_boundary_.insert(std::make_pair(5, left_boundary_.at(7) - temp + (vx / 6.0 + vy / 2.0) * rho));
		left_boundary_.insert(std::make_pair(8, left_boundary_.at(6) + temp + (vx / 6.0 - vy / 2.0) * rho));

		// + set velocities to fluid->size().second - 1 row after recalculation
		for (auto id : ids_1)
			left_boundary_.erase(id);

		for (auto id : ids_2)
			left_boundary_.erase(id);

	}
	else if (first == Boundary::RIGHT)
	{
		std::vector<double> rho(y_size, 0.0);

		for (auto id : ids_1)
			rho = rho + right_boundary_.at(id);
		for (auto id : ids_2)
			rho = rho + right_boundary_.at(id) * 2.0;

		rho = rho / (1.0 + vx);

		right_boundary_.insert(std::make_pair(3, right_boundary_.at(1) - 2.0 / 3.0 * rho * vx));


		std::vector<double> temp = (right_boundary_.at(4) - right_boundary_.at(2)) / 2.0;

		right_boundary_.insert(std::make_pair(6, right_boundary_.at(8) + temp - (vx / 6.0 - vy / 2.0) * rho));
		right_boundary_.insert(std::make_pair(7, right_boundary_.at(5) - temp + (vx / 6.0 + vy / 2.0) * rho));

		// + set velocities to fluid->size().second - 1 row after recalculation
		for (auto id : ids_1)
			right_boundary_.erase(id);

		for (auto id : ids_2)
			right_boundary_.erase(id);

	}

}

void BCs::SwapId(std::map<int, std::vector<double>> & map, int const from, int const to)
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

#pragma region 3d

#include"../../math/array_func_impl.h"

BCs3D::BCs3D(int rows, int colls, DistributionFunction3D<double>& dfunc) :
	height_(rows), length_(colls - 2), f_ptr_(&dfunc)
{
	
}

bool BCs3D::PrepareValuesForSingleBC(Boundary const BC, BCType const bc_type)
{
	// Pointer to function, which gets appropriate values, bepending on wall type:
	// for TOP : ptr to GetTopBoundaryValues, for BOTTOM ptr to GetBottomBoundaryValue
	std::vector<double>(DistributionFunction3D<double>::*ptrToFunc)(int) const = nullptr;

	switch (BC)
	{
	case Boundary::TOP:

		// !!! VON NEUMANN IMPLEMENTATION PROCESS
		ptrToFunc = &DistributionFunction3D<double>::GetTopBoundaryValues;

		if(bc_type == BCType::BOUNCE_BACK || bc_type == BCType::PERIODIC)
			return WriteBoundaryValues(bc_type, top_boundary_, top_ids_, ptrToFunc);
		else if (bc_type == BCType::VON_NEUMAN)
			return WriteVonNeumannBoundaryValues(bc_type, top_boundary_, middle_layer_ids_, top_ids_, ptrToFunc);
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
	/* Pointer to function to distinguish GETTING boundaries: TOP, BOTTOM, e.t.c */ std::vector<double> (DistributionFunction3D<double>::*ptrToFunc)(int) const)
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

bool BCs3D::WriteVonNeumannBoundaryValues(BCType const bc_type, std::map<int, std::vector<double>>& bc_boundary, const std::vector<int>& bc_ids_1, const std::vector<int>& bc_ids_2, std::vector<double>(DistributionFunction3D<double>::* ptrToFunc)(int) const)
{
	if (bc_type == BCType::VON_NEUMAN)
	{
		// Idea: Get Values from first and second layers: for example for BOTTOM we need ALL distribution function from MIDDLE and TOP layers
		// Get appropriate values from probability distribution function (using APPROPRIATE function via ptr. to function) 
		// and write them in appropriate class fields.

		for (auto bc_id_1 : bc_ids_1)
			bc_boundary.insert(std::make_pair(bc_id_1, (*f_ptr_.*ptrToFunc)(bc_id_1)));
		for (auto bc_id_2 : bc_ids_2)
			bc_boundary.insert(std::make_pair(bc_id_2, (*f_ptr_.*ptrToFunc)(bc_id_2)));

		return true;
	}
	else
	{

	}
	return false;
}

void BCs3D::PrepareValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc, BCType const near_bc, BCType far_bc)
{
	if (PrepareValuesForSingleBC(Boundary::TOP, top_bc) &&
		PrepareValuesForSingleBC(Boundary::BOTTOM, bottm_bc) &&
		PrepareValuesForSingleBC(Boundary::LEFT, left_bc) &&
		PrepareValuesForSingleBC(Boundary::RIGHT, right_bc) && 
		PrepareValuesForSingleBC(Boundary::NEAR, near_bc) &&
		PrepareValuesForSingleBC(Boundary::FAAR, far_bc) )
	{
		// Лог что все значения получилось взять
	}
	else
	{
		std::cout << "Could not prepare values for all boundaries.\n";
		throw;
	}
}

bool BCs3D::RecordValuesForSingleBC(Boundary const BC, BCType const bc_type)
{
	void(DistributionFunction3D<double>::*ptrToFunc)(int const, std::vector<double> const &) = nullptr;

	switch (BC)
	{
	case Boundary::TOP:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetBottomBoundaryValue;
			return RecordBoundaryValues(bc_type, bottom_boundary_, top_ids_, ptrToFunc);
		}
		else if (bc_type == BCType::BOUNCE_BACK)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetTopBoundaryValue;
			return RecordBoundaryValues(bc_type, top_boundary_, bottom_ids_, ptrToFunc);
		}
		else if (bc_type == BCType::VON_NEUMAN)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetTopBoundaryValue;
			return RecordBoundaryValues(bc_type, top_boundary_, bottom_ids_, ptrToFunc);
		}
		break;

	case Boundary::BOTTOM:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetTopBoundaryValue;
			return RecordBoundaryValues(bc_type, top_boundary_, bottom_ids_, ptrToFunc);

		}
		else if (bc_type == BCType::BOUNCE_BACK)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetBottomBoundaryValue;
			return RecordBoundaryValues(bc_type, bottom_boundary_, top_ids_, ptrToFunc);
		}
		break;

	case Boundary::RIGHT:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetLeftBoundaryValue;
			return RecordBoundaryValues(bc_type, left_boundary_, right_ids_, ptrToFunc);
		}
		else if (bc_type == BCType::BOUNCE_BACK)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetRightBoundaryValue;
			return RecordBoundaryValues(bc_type, right_boundary_, left_ids_, ptrToFunc);
		}

		break;

	case Boundary::LEFT:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetRightBoundaryValue;
			return RecordBoundaryValues(bc_type, right_boundary_, left_ids_, ptrToFunc);
		}
		else if (bc_type == BCType::BOUNCE_BACK)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetLeftBoundaryValue;
			return RecordBoundaryValues(bc_type, left_boundary_, right_ids_, ptrToFunc);
		}
		break;

	case Boundary::NEAR:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetFarBoundaryValue;
			return RecordBoundaryValues(bc_type, far_boundary_, near_ids_, ptrToFunc);

		}
		else if (bc_type == BCType::BOUNCE_BACK)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetNearBoundaryValue;
			return RecordBoundaryValues(bc_type, near_boundary_, far_ids_, ptrToFunc);
		}
		break;

	case Boundary::FAAR:

		if (bc_type == BCType::PERIODIC)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetNearBoundaryValue;
			return RecordBoundaryValues(bc_type, near_boundary_, far_ids_, ptrToFunc);
		}
		else if (bc_type == BCType::BOUNCE_BACK)
		{
			ptrToFunc = &DistributionFunction3D<double>::SetFarBoundaryValue;
			return RecordBoundaryValues(bc_type, far_boundary_, near_ids_, ptrToFunc);
		}
		break;

	default:
		std::cout << "Wrong Boundary type is used while prepair values for BC.\n";
		return false;
		break;
	}
}

bool BCs3D::RecordBoundaryValues(BCType const bc_type, std::map<int, std::vector<double>>& bc_boundary, const std::vector<int>& bc_ids, 
	/* Pointer to function to distinguish GETTING boundaries: TOP, BOTTOM, e.t.c */ void(DistributionFunction3D<double>::* ptrToFunc)(int const, std::vector<double> const &))
{
	if (bc_type == BCType::PERIODIC || bc_type == BCType::BOUNCE_BACK)
	{
		for (auto bc_id : bc_ids)
		{
			auto f_id = bc_boundary.find(bc_id);
			(*f_ptr_.*ptrToFunc)(f_id->first, f_id->second);
		}
		return true;
	}
	else if (bc_type == BCType::VON_NEUMAN)
	{
		// !!! In testing mode
		for (auto bc_id : bc_ids)
		{
			auto f_id = bc_boundary.find(bc_id);
			(*f_ptr_.*ptrToFunc)(f_id->first, f_id->second);
		}
		return true;
	}
	else
	{
		std::cout << "Wrong BC Type is used!\n";
	}
	return false;
}


void BCs3D::RecordValuesForAllBC(BCType const top_bc, BCType const bottm_bc, BCType const left_bc, BCType const right_bc, BCType const near_bc, BCType far_bc)
{
	if (RecordValuesForSingleBC(Boundary::TOP, top_bc) &&
		RecordValuesForSingleBC(Boundary::BOTTOM, bottm_bc) &&
		RecordValuesForSingleBC(Boundary::LEFT, left_bc) &&
		RecordValuesForSingleBC(Boundary::RIGHT, right_bc) &&
		RecordValuesForSingleBC(Boundary::NEAR, near_bc) &&
		RecordValuesForSingleBC(Boundary::FAAR, far_bc))
	{
		f_ptr_->ClearBoundaries();
	}
	else
	{
		std::cout << "Could not record values for all boundaries.\n";
		throw;
	}

	
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
		SwapIds(top_boundary_, 9, 14);
		SwapIds(top_boundary_, 10, 17);
		SwapIds(top_boundary_, 11, 18);
		SwapIds(top_boundary_, 12, 15);
		SwapIds(top_boundary_, 13, 16);
	}
	else if (first == Boundary::BOTTOM) 
	{
		SwapIds(bottom_boundary_, 14, 9);
		SwapIds(bottom_boundary_, 15, 12);
		SwapIds(bottom_boundary_, 16, 13);
		SwapIds(bottom_boundary_, 17, 10);
		SwapIds(bottom_boundary_, 18, 11);
	}
	else if (first == Boundary::LEFT) 
	{
		SwapIds(left_boundary_, 3, 1);
		SwapIds(left_boundary_, 6, 8);
		SwapIds(left_boundary_, 7, 5);
		SwapIds(left_boundary_, 12, 15);
		SwapIds(left_boundary_, 17, 10);
	}
	else if (first == Boundary::RIGHT) 
	{
		SwapIds(right_boundary_, 1, 3);
		SwapIds(right_boundary_, 5, 7);
		SwapIds(right_boundary_, 8, 6);
		SwapIds(right_boundary_, 10, 17);
		SwapIds(right_boundary_, 15, 12);
	}
	else if (first == Boundary::NEAR)
	{
		SwapIds(near_boundary_, 4, 2);
		SwapIds(near_boundary_, 7, 5);
		SwapIds(near_boundary_, 8, 6);
		SwapIds(near_boundary_, 13, 16);
		SwapIds(near_boundary_, 18, 11);
	}
	else if (first == Boundary::FAAR)
	{
		SwapIds(far_boundary_, 2, 4);
		SwapIds(far_boundary_, 5, 7);
		SwapIds(far_boundary_, 6, 8);
		SwapIds(far_boundary_, 11, 18);
		SwapIds(far_boundary_, 16, 13);
	}
}

void BCs3D::VonNeumannBC(Boundary const first, const double vx, const double vy, const double vz)
{
	if (first == Boundary::TOP)
	{
		typedef std::vector< std::pair<int, std::vector<double>>> VectorOfPairs;

		// Copy from map to vector for easy work
		const int size = top_boundary_.begin()->second.size();
		std::vector<double> rho(size, 0.0);
		// >>> Calculate rho
		for (auto middleId : middle_layer_ids_)
			rho = rho + top_boundary_.at(middleId);

		for (auto topId : top_ids_)
				rho = rho + top_boundary_.at(topId) * 2.0;
			
		rho = rho / (1.0 + vz);
		// >>>
		

		// >>> Calculate  coefs N
		std::vector<double> Nxz(size, 0.0);
		Nxz = 0.5 * (top_boundary_.at(1) + top_boundary_.at(5) + top_boundary_.at(8) - (top_boundary_.at(3) + top_boundary_.at(6) + top_boundary_.at(7))) - 1.0 / 3.0 * rho * vx;

		std::vector<double> Nyz(size, 0.0);
		Nyz = 0.5 * (top_boundary_.at(2) + top_boundary_.at(5) + top_boundary_.at(6) - (top_boundary_.at(4) + top_boundary_.at(7) + top_boundary_.at(8))) - 1.0 / 3.0 * rho * vy;

		// >>>

		auto f14 = top_boundary_.at(9) - 1.0 / 3.0 * rho * vz;
		auto f15 = top_boundary_.at(12) + rho / 6.0 *(-vz + vx) - Nxz;
		auto f17 = top_boundary_.at(10) + rho / 6.0 *(-vz - vx) + Nxz;
		auto f16 = top_boundary_.at(13) + rho / 6.0 *(-vz + vx) - Nyz; // May be change signs
		auto f18 = top_boundary_.at(11) + rho / 6.0 *(-vz - vx) + Nyz; // May be change signs

		top_boundary_.clear();
		top_boundary_.insert(std::make_pair(14, f14));
		top_boundary_.insert(std::make_pair(15, f15));
		top_boundary_.insert(std::make_pair(16, f16));
		top_boundary_.insert(std::make_pair(17, f17));
		top_boundary_.insert(std::make_pair(18, f18));
	}
	else
	{
		std::cout << "Is not implemented yet!";
		throw;
	}
}

void BCs3D::SwapIds(std::map<int, std::vector<double>>& map, int const from, int const to)
{
	std::vector<double> temp;
	auto iter = map.find(from);
	temp.swap(iter->second);
	map.erase(iter);
	map.insert(std::make_pair(to, temp));
}

#pragma endregion