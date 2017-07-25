#pragma once

#include<memory>
#include<map>

#include"../../modeling_area/fluid.h"
#include"../../math/array_func_impl.h" // Для работы с ГУ Фон-Неймана

class MyBCs
{
public:
	MyBCs(Boundary bound, Fluid* fluid, Medium & medium) : bound_(bound), fluid_(fluid), medium_(medium), ptrToGetFunc(nullptr), ptrToSetFunc(nullptr) {}
	//virtual ~MyBCs() = 0;

	virtual void ApplyBCs() {}

protected:
	Fluid* fluid_;
	Medium& medium_;
	std::map<int, std::vector<double>> dfunc_;
	Boundary bound_;

	std::vector<double> (DistributionFunction<double>::* ptrToGetFunc)(int const) const;
	void(DistributionFunction<double>::* ptrToSetFunc)(int const, std::vector<double> const &);

	std::vector<int> prepare_ids_;
	std::vector<int> record_ids_;

	const std::map<Boundary, std::vector<int>> ids_
	{
		{ Boundary::TOP, { 2, 5, 6 } },
		{ Boundary::BOTTOM, { 4, 7, 8 } },
		{ Boundary::LEFT, { 3, 6, 7 } },
		{ Boundary::RIGHT, { 1, 8, 5 } },
	};

	const std::vector<int> mid_h_ids_{ 0, 1, 3 };
	const std::vector<int> mid_v_ids_{ 0, 2, 4 };
};


class PeriodicBCs : public MyBCs
{
public:
	PeriodicBCs(Boundary bound, Fluid* fluid, Medium & medium) : MyBCs(bound, fluid, medium)
	{
		prepare_ids_ = ids_.at(bound);
		record_ids_ = ids_.at(bound);

		switch (bound)
		{
		case Boundary::TOP:
			ptrToGetFunc = &DistributionFunction<double>::GetTopValues;
			ptrToSetFunc = &DistributionFunction<double>::SetBottomValues;
			break;
		case Boundary::BOTTOM:
			ptrToGetFunc = &DistributionFunction<double>::GetBottomValues;
			ptrToSetFunc = &DistributionFunction<double>::SetTopValues;
			break;
		case Boundary::RIGHT:
			ptrToGetFunc = &DistributionFunction<double>::GetRightValues;
			ptrToSetFunc = &DistributionFunction<double>::SetLeftValues;
			break;
		case Boundary::LEFT:
			ptrToGetFunc = &DistributionFunction<double>::GetLeftValues;
			ptrToSetFunc = &DistributionFunction<double>::SetRightValues;
			break;
		default:
			std::cout << "Wrong input BC boundary type for Periodic BC.\n";
			break;
		}
	}

	~PeriodicBCs() {}

	void ApplyBCs() override
	{
		for (int q = 0; q < 3; ++q)
		{
			std::vector<double> temp = (fluid_->f_.*ptrToGetFunc)(prepare_ids_.at(q));
			(fluid_->f_.*ptrToSetFunc)(record_ids_.at(q), temp);
		}
	}
};


class BounceBackBCs: public MyBCs
{
public:
	BounceBackBCs(Boundary bound, Fluid* fluid, Medium & medium) : MyBCs(bound, fluid, medium)
	{
		prepare_ids_ = ids_.at(bound);

		switch (bound)
		{
		case Boundary::TOP:
			ptrToGetFunc = &DistributionFunction<double>::GetTopValues;
			ptrToSetFunc = &DistributionFunction<double>::SetTopValues;

			record_ids_ = ids_.at(Boundary::BOTTOM);
			break;
		case Boundary::BOTTOM:
			ptrToGetFunc = &DistributionFunction<double>::GetBottomValues;
			ptrToSetFunc = &DistributionFunction<double>::SetBottomValues;

			record_ids_ = ids_.at(Boundary::TOP);
			break;
		case Boundary::RIGHT:
			ptrToGetFunc = &DistributionFunction<double>::GetRightValues;
			ptrToSetFunc = &DistributionFunction<double>::SetRightValues;

			record_ids_ = ids_.at(Boundary::LEFT);
			break;
		case Boundary::LEFT:
			ptrToGetFunc = &DistributionFunction<double>::GetLeftValues;
			ptrToSetFunc = &DistributionFunction<double>::SetLeftValues;

			record_ids_ = ids_.at(Boundary::RIGHT);
			break;
		default:
			std::cout << "Wrong input BC boundary type for Periodic BC.\n";
			break;
		}
	}

	~BounceBackBCs() {}

	void ApplyBCs() override
	{
		for (int q = 0; q < 3; ++q)
		{
			std::vector<double> temp = (fluid_->f_.*ptrToGetFunc)(prepare_ids_.at(q));
			(fluid_->f_.*ptrToSetFunc)(record_ids_.at(q), temp);
		}
	}
};

class VonNeumannBCs : public MyBCs
{
public:

	VonNeumannBCs(Boundary bound, Fluid* fluid, Medium & medium, std::vector<double> & vx, std::vector<double> & vy) : MyBCs(bound, fluid, medium), vx_(vx), vy_(vy)
	{
		const int rows = fluid_->GetRowsNumber();
		const int colls = fluid->GetColumnsNumber();

		// Make assert on size of arrays of vx, vy
		prepare_ids_ = ids_.at(bound);

		switch (bound)
		{
		case Boundary::TOP:
			ptrToGetFunc = &DistributionFunction<double>::GetTopValues1;
			ptrToSetFunc = &DistributionFunction<double>::SetTopValues;

			record_ids_ = ids_.at(Boundary::BOTTOM);
			mid_ids_ = mid_h_ids_;
			break;
		case Boundary::BOTTOM:
			ptrToGetFunc = &DistributionFunction<double>::GetBottomValues1;
			ptrToSetFunc = &DistributionFunction<double>::SetBottomValues;

			record_ids_ = ids_.at(Boundary::TOP);
			mid_ids_ = mid_h_ids_;
			break;
		case Boundary::RIGHT:
			ptrToGetFunc = &DistributionFunction<double>::GetRightValues1;
			ptrToSetFunc = &DistributionFunction<double>::SetRightValues;

			record_ids_ = ids_.at(Boundary::LEFT);
			mid_ids_ = mid_v_ids_;
			break;
		case Boundary::LEFT:
			ptrToGetFunc = &DistributionFunction<double>::GetLeftValues1;
			ptrToSetFunc = &DistributionFunction<double>::SetLeftValues;

			record_ids_ = ids_.at(Boundary::RIGHT);
			mid_ids_ = mid_v_ids_;
			break;
		default:
			std::cout << "Wrong input BC boundary type for Periodic BC.\n";
			break;
		}


		for (int q = 0; q < kQ; ++q)
		{
			corner_rt_.at(q) = fluid->f_.Get(q, rows - 2, colls - 2);
			corner_rb_.at(q) = fluid->f_.Get(q, rows - 2, 1);
			corner_lt_.at(q) = fluid->f_.Get(q, 1, colls - 2);
			corner_lb_.at(q) = fluid->f_.Get(q, rows - 2, 1);
		}
	}
	~VonNeumannBCs() {}

	void ApplyBCs() override
	{
		const int size = (bound_ == Boundary::TOP || bound_ == Boundary::BOTTOM) ? fluid_->GetColumnsNumber() - 2 : fluid_->GetRowsNumber() - 2;
		std::vector<double> rho(size, 0.0);
		std::vector<double> temp;
		std::vector<double> cur_f;

		for (const auto & q : mid_ids_)
			rho = rho + (fluid_->f_.*ptrToGetFunc)(q);
		for (const auto & q : prepare_ids_)
			rho = rho +  (fluid_->f_.*ptrToGetFunc)(q) * 2.0;

		switch (bound_)
		{
		case Boundary::TOP:
			rho = rho / (1.0 + vy_);

			cur_f = (fluid_->f_.*ptrToGetFunc)(2) - 2.0 / 3.0 * rho * vy_;
			(fluid_->f_.*ptrToSetFunc)(4, cur_f);

			temp = ((fluid_->f_.*ptrToGetFunc)(1) - (fluid_->f_.*ptrToGetFunc)(3)) / 2.0;
			//f7
			cur_f = (fluid_->f_.*ptrToGetFunc)(5) + temp - (vy_ / 6.0 + vx_ / 2.0) * rho;
			(fluid_->f_.*ptrToSetFunc)(7, cur_f);
			//f8
			cur_f = (fluid_->f_.*ptrToGetFunc)(6) - temp - (vy_ / 6.0 - vx_ / 2.0) * rho;
			(fluid_->f_.*ptrToSetFunc)(8, cur_f);
			break;
		case Boundary::BOTTOM:
			rho = rho / (1.0 - vy_);
			// f2
			cur_f = (fluid_->f_.*ptrToGetFunc)(4) + 2.0 / 3.0 * rho * vy_;
			(fluid_->f_.*ptrToSetFunc)(2, cur_f);

			temp = ((fluid_->f_.*ptrToGetFunc)(1) - (fluid_->f_.*ptrToGetFunc)(3)) / 2.0;
			//f5
			cur_f = (fluid_->f_.*ptrToGetFunc)(7) - temp + (vy_ / 6.0 + vx_ / 2.0) * rho;
			(fluid_->f_.*ptrToSetFunc)(5, cur_f);
			//f6
			cur_f = (fluid_->f_.*ptrToGetFunc)(8) + temp + (vy_ / 6.0 - vx_ / 2.0) * rho;
			(fluid_->f_.*ptrToSetFunc)(6, cur_f);

			break;
		case Boundary::RIGHT:
			rho = rho / (1.0 + vx_);
			// f1
			cur_f = (fluid_->f_.*ptrToGetFunc)(1) - 2.0 / 3.0 * rho * vx_;
			(fluid_->f_.*ptrToSetFunc)(3, cur_f);

			temp = ((fluid_->f_.*ptrToGetFunc)(2) - (fluid_->f_.*ptrToGetFunc)(4)) / 2.0;
			// f5
			cur_f = (fluid_->f_.*ptrToGetFunc)(8) - temp - (vx_ / 6.0 - vy_ / 2.0) * rho;
			(fluid_->f_.*ptrToSetFunc)(6, cur_f);
			// f6
			cur_f = (fluid_->f_.*ptrToGetFunc)(5) + temp - (vx_ / 6.0 + vy_ / 2.0) * rho;
			(fluid_->f_.*ptrToSetFunc)(7, cur_f);

			break;
		case Boundary::LEFT:

			rho = rho / (1.0 - vx_);
			// f1
			cur_f = (fluid_->f_.*ptrToGetFunc)(3) + 2.0 / 3.0 * rho * vx_;
			(fluid_->f_.*ptrToSetFunc)(1, cur_f);

			temp = ((fluid_->f_.*ptrToGetFunc)(2) - (fluid_->f_.*ptrToGetFunc)(4)) / 2.0;
			// f5
			cur_f = (fluid_->f_.*ptrToGetFunc)(7) - temp + (vx_ / 6.0 + vy_ / 2.0) * rho;
			(fluid_->f_.*ptrToSetFunc)(5, cur_f);
			// f6
			cur_f = (fluid_->f_.*ptrToGetFunc)(6) + temp + (vx_ / 6.0 - vy_ / 2.0) * rho;
			(fluid_->f_.*ptrToSetFunc)(8, cur_f);
			break;
		default:
			break;
		}

	
	}

private:
	std::vector<int> mid_ids_;
	std::vector<double> vx_;
	std::vector<double> vy_;


	std::array<double, kQ> corner_rt_;
	std::array<double, kQ> corner_rb_;
	std::array<double, kQ> corner_lt_;
	std::array<double, kQ> corner_lb_;
};

class DirichletBCs : public MyBCs
{
public:
	DirichletBCs(Boundary bound, Fluid* fluid, Medium & medium, std::vector<double> & rho) : MyBCs(bound, fluid, medium), rho_(rho)
	{
		const int rows = fluid_->GetRowsNumber();
		const int colls = fluid->GetColumnsNumber();

		prepare_ids_ = ids_.at(bound);

		switch (bound)
		{
		case Boundary::TOP:
			ptrToGetFunc = &DistributionFunction<double>::GetTopValues1;
			ptrToSetFunc = &DistributionFunction<double>::SetTopValues;

			record_ids_ = ids_.at(Boundary::BOTTOM);
			mid_ids_ = mid_h_ids_;
			break;
		case Boundary::BOTTOM:
			ptrToGetFunc = &DistributionFunction<double>::GetBottomValues1;
			ptrToSetFunc = &DistributionFunction<double>::SetBottomValues;

			record_ids_ = ids_.at(Boundary::TOP);
			mid_ids_ = mid_h_ids_;
			break;
		case Boundary::RIGHT:
			ptrToGetFunc = &DistributionFunction<double>::GetRightValues1;
			ptrToSetFunc = &DistributionFunction<double>::SetRightValues;

			record_ids_ = ids_.at(Boundary::LEFT);
			mid_ids_ = mid_v_ids_;
			break;
		case Boundary::LEFT:
			ptrToGetFunc = &DistributionFunction<double>::GetLeftValues1;
			ptrToSetFunc = &DistributionFunction<double>::SetLeftValues;

			record_ids_ = ids_.at(Boundary::RIGHT);
			mid_ids_ = mid_v_ids_;
			break;
		default:
			std::cout << "Wrong input BC boundary type for Periodic BC.\n";
			break;
		}
	}
	~DirichletBCs() {}

	void ApplyBCs() override
	{
		const int size = (bound_ == Boundary::TOP || bound_ == Boundary::BOTTOM) ? fluid_->GetColumnsNumber() - 2 : fluid_->GetRowsNumber() - 2;
		std::vector<double> v(size, 0.0);
		std::vector<double> temp;
		std::vector<double> cur_f;

		for (const auto & q : mid_ids_)
			v = v + (fluid_->f_.*ptrToGetFunc)(q);
		for (const auto & q : prepare_ids_)
			v = v + (fluid_->f_.*ptrToGetFunc)(q) * 2.0;

		switch (bound_)
		{
		case Boundary::TOP:
			v = v / rho_ - 1.0;

			cur_f = (fluid_->f_.*ptrToGetFunc)(2) - 2.0 / 3.0 * rho_ * v;
			(fluid_->f_.*ptrToSetFunc)(4, cur_f);

			temp = ((fluid_->f_.*ptrToGetFunc)(1) - (fluid_->f_.*ptrToGetFunc)(3)) / 2.0;
			//f7
			cur_f = (fluid_->f_.*ptrToGetFunc)(5) + temp - v * rho_ / 6.0;
			(fluid_->f_.*ptrToSetFunc)(7, cur_f);
			//f8
			cur_f = (fluid_->f_.*ptrToGetFunc)(6) - temp - v * rho_ / 6.0;
			(fluid_->f_.*ptrToSetFunc)(8, cur_f);
			break;
		case Boundary::BOTTOM:
			v = (v / rho_) * -1.0 + 1.0;
			// f2
			cur_f = (fluid_->f_.*ptrToGetFunc)(4) + 2.0 / 3.0 * rho_ * v;
			(fluid_->f_.*ptrToSetFunc)(2, cur_f);

			temp = ((fluid_->f_.*ptrToGetFunc)(1) - (fluid_->f_.*ptrToGetFunc)(3)) / 2.0;
			//f5
			cur_f = (fluid_->f_.*ptrToGetFunc)(7) - temp + rho_ * v / 6.0;
			(fluid_->f_.*ptrToSetFunc)(5, cur_f);
			//f6
			cur_f = (fluid_->f_.*ptrToGetFunc)(8) + temp + rho_ * v / 6.0;
			(fluid_->f_.*ptrToSetFunc)(6, cur_f);

			break;
		case Boundary::RIGHT:
			v = v / rho_ - 1.0;
			// f1
			cur_f = (fluid_->f_.*ptrToGetFunc)(1) - 2.0 / 3.0 * rho_ * v;
			(fluid_->f_.*ptrToSetFunc)(3, cur_f);

			temp = ((fluid_->f_.*ptrToGetFunc)(2) - (fluid_->f_.*ptrToGetFunc)(4)) / 2.0;
			// f5
			cur_f = (fluid_->f_.*ptrToGetFunc)(8) - temp - rho_ * v / 6.0 ;
			(fluid_->f_.*ptrToSetFunc)(6, cur_f);
			// f6
			cur_f = (fluid_->f_.*ptrToGetFunc)(5) + temp - rho_ * v / 6.0;
			(fluid_->f_.*ptrToSetFunc)(7, cur_f);

			break;
		case Boundary::LEFT:
			v = (v / rho_) * -1.0 + 1.0;
			// f1
			cur_f = (fluid_->f_.*ptrToGetFunc)(3) + 2.0 / 3.0 *rho_ * v;
			(fluid_->f_.*ptrToSetFunc)(1, cur_f);

			temp = ((fluid_->f_.*ptrToGetFunc)(2) - (fluid_->f_.*ptrToGetFunc)(4)) / 2.0;
			// f5
			cur_f = (fluid_->f_.*ptrToGetFunc)(7) - temp + rho_ * v / 6.0;
			(fluid_->f_.*ptrToSetFunc)(5, cur_f);
			// f6
			cur_f = (fluid_->f_.*ptrToGetFunc)(6) + temp + rho_ * v / 6.0;
			(fluid_->f_.*ptrToSetFunc)(8, cur_f);
			break;
		default:
			break;
		}

	}


private:
	std::vector<int> mid_ids_;
	std::vector<double> rho_;
	
};

