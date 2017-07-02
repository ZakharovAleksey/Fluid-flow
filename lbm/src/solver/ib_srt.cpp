#include"ib_srt.h"

IBSolver::IBSolver(double tau, Fluid & fluid, Medium & medium, std::unique_ptr<ImmersedBody> body) : tau_(tau)
{
	CreateDataFolder("Data");
	CreateDataFolder("Data\\ib_lbm_data");
	CreateDataFolder("Data\\ib_lbm_data\\body_form_txt");
	CreateDataFolder("Data\\ib_lbm_data\\body_form_vtk");
	CreateDataFolder("Data\\ib_lbm_data\\fluid_txt");
	CreateDataFolder("Data\\ib_lbm_data\\fluid_vtk");


	std::cout << " --- Input parameters :\n";
	std::cout << "nu = " << (tau - 0.5) / 3.0 << std::endl;
	std::cout << "Re =  " << fluid_->size().second * 0.001 / ((tau - 0.5) / 3.0) << std::endl;

	fluid_ = std::unique_ptr<Fluid>(new Fluid(fluid));
	medium_ = std::unique_ptr<Medium>(new Medium(medium));
	body_ = std::move(body);//std::unique_ptr<ImmersedBody>(new ImmersedBody(body));

	int rows = fluid_->size().first;
	int colls = fluid_->size().second;

	fx_ = std::make_unique<Matrix2D<double>>(rows, colls);
	fy_ = std::make_unique<Matrix2D<double>>(rows, colls);

	force_member_.resize(kQ, 0.0);
}

IBSolver::IBSolver(double tau, Fluid & fluid, Medium & medium, std::vector<ImmersedBody*> bodies) : tau_(tau)
{
	CreateDataFolder("Data");
	CreateDataFolder("Data\\ib_lbm_data");
	CreateDataFolder("Data\\ib_lbm_data\\body_form_txt");
	CreateDataFolder("Data\\ib_lbm_data\\body_form_vtk");
	CreateDataFolder("Data\\ib_lbm_data\\fluid_txt");
	CreateDataFolder("Data\\ib_lbm_data\\fluid_vtk");

	std::cout << " --- Input parameters :\n";
	std::cout << "nu = " << (tau - 0.5) / 3.0 << std::endl;
	std::cout << "Re =  " << fluid.size().second * 0.01 / ((tau - 0.5) / 3.0) << std::endl;
	fluid_ = std::unique_ptr<Fluid>(new Fluid(fluid));
	medium_ = std::unique_ptr<Medium>(new Medium(medium));

	int i = 0;
	for (const auto& e : bodies)
		im_bodies_.push_back(bodies.at(i++));

	int rows = fluid_->size().first;
	int colls = fluid_->size().second;

	fx_ = std::make_unique<Matrix2D<double>>(rows, colls);
	fy_ = std::make_unique<Matrix2D<double>>(rows, colls);

	force_member_.resize(kQ, 0.0);
}

void IBSolver::feqCalculate()
{
	// ѕроверить надо ли, или без нее все нормально
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

void IBSolver::Streaming()
{
	for (int q = 0; q < kQ; ++q)
	{
		Matrix2D<double> temp = fluid_->f_[q];
		fluid_->f_[q].FillWith(0.0);

		for (unsigned y = 0; y < fluid_->size().first; ++y)
			for (unsigned x = 0; x < fluid_->size().second; ++x)
				if (medium_->IsFluid(y, x))
					fluid_->f_[q](y + kEy[q], x + kEx[q]) = temp(y, x);
	}

	fluid_->f_.fillBoundaries(0.0);
}

void IBSolver::CalculateForces()
{
	double gravity = 0.0;

	for (int y = 0; y < fluid_->size().first; ++y)
		for (int x = 0; x < fluid_->size().second; ++x)
		{
			force_member_.at(0) = (1.0 - 0.5 / tau_) * kW[0] * (3.0 * ((-fluid_->vx_(y, x) * ((*fx_)(y, x) + gravity) + -fluid_->vy_(y, x) * (*fy_)(y, x))));
			force_member_[1] = (1 - 0.5 / tau_) * kW[1] * (3.0 * ((1 - fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity) + (-fluid_->vy_(y, x)) * (*fy_)(y, x)) + 9.0 * (fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity));
			force_member_[2] = (1 - 0.5 / tau_) * kW[2] * (3.0* ((-1 - fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity) + (-fluid_->vy_(y, x)) * (*fy_)(y, x)) + 9.0 * (fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity));
			force_member_[3] = (1 - 0.5 / tau_) * kW[3] * (3.0 * ((-fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity) + (1 - fluid_->vy_(y, x)) * (*fy_)(y, x)) + 9.0 * (fluid_->vy_(y, x)) * (*fy_)(y, x));
			force_member_[4] = (1 - 0.5 / tau_) * kW[4] * (3.0 * ((-fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity) + (-1 - fluid_->vy_(y, x)) * (*fy_)(y, x)) + 9.0 * (fluid_->vy_(y, x)) * (*fy_)(y, x));
			force_member_[5] = (1 - 0.5 / tau_) * kW[5] * (3.0 * ((1 - fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity) + (1 - fluid_->vy_(y, x)) * (*fy_)(y, x)) + 9.0 * (fluid_->vx_(y, x) + fluid_->vy_(y, x)) * ((*fx_)(y, x) + gravity + (*fy_)(y, x)));
			force_member_[6] = (1 - 0.5 / tau_) * kW[6] * (3.0 * ((-1 - fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity) + (-1 - fluid_->vy_(y, x)) * (*fy_)(y, x)) + 9.0 * (fluid_->vx_(y, x) + fluid_->vy_(y, x)) * ((*fx_)(y, x) + gravity + (*fy_)(y, x)));
			force_member_[7] = (1 - 0.5 / tau_) * kW[7] * (3.0 * ((1 - fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity) + (-1 - fluid_->vy_(y, x)) * (*fy_)(y, x)) + 9.0 * (fluid_->vx_(y, x) - fluid_->vy_(y, x)) * ((*fx_)(y, x) + gravity - (*fy_)(y, x)));
			force_member_[8] = (1 - 0.5 / tau_) * kW[8] * (3.0 * ((-1 - fluid_->vx_(y, x)) * ((*fx_)(y, x) + gravity) + (1 - fluid_->vy_(y, x)) * (*fy_)(y, x)) + 9 * (fluid_->vx_(y, x) - fluid_->vy_(y, x)) * ((*fx_)(y, x) + gravity - (*fy_)(y, x)));
		}
}

void IBSolver::Collision()
{
	CalculateForces();

	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] += (fluid_->feq_[q] - fluid_->f_[q]) / tau_ + force_member_.at(q);
}

void IBSolver::Recalculate()
{
	fluid_->rho_ = fluid_->f_.calculateDensity();

	// Recalculate velocities with additional term
	fluid_->vx_ = fluid_->f_.calculateVelocity(kEx, fluid_->rho_, *fx_);
	fluid_->vy_ = fluid_->f_.calculateVelocity(kEy, fluid_->rho_, *fy_);
}

void IBSolver::Solve(int iter_numb)
{
	feqCalculate();

	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] = fluid_->feq_[q];

	BCs BC(fluid_->f_);
	Microphone mic;

	for (int iter = 0; iter < iter_numb; ++iter)
	{
		//mic.PerformMeasurements(iter, im_bodies_, fluid_->vx_, "vx");
		//mic.PerformMeasurements(iter, im_bodies_, fluid_->vx_, "vy");
		//mic.PerformMeasurements(iter, im_bodies_, fluid_->vx_, "rho");

		// Clean fx, fy fields
		fx_->FillWith(0.0);
		fy_->FillWith(0.0);

		for (auto& i : im_bodies_)
		{
			i->CalculateForces();
			//i->SpreadForces(*fx_, *fy_);
		}

		// Deal with RBC-Wall intearaction
		//Interaction(im_bodies_.at(2), im_bodies_.at(0));
		//Interaction(im_bodies_.at(2), im_bodies_.at(1));

		for (auto& i : im_bodies_)
		{
			//i->CalculateForces();
			i->SpreadForces(*fx_, *fy_);
		}

		Collision();
		BC.PrepareValuesForAllBC(BCType::PERIODIC, BCType::PERIODIC, BCType::VON_NEUMAN, BCType::VON_NEUMAN);
		if (medium_->IsImmersedBodies())
			BC.PrepareAdditionalBCs(*medium_);
		

		Streaming();

		//BC.PrepareAdditionalBCs(*medium_);

		BC.PeriodicBC(Boundary::TOP, Boundary::BOTTOM);
		//BC.BounceBackBC(Boundary::TOP);
		//BC.BounceBackBC(Boundary::BOTTOM);

		double vx = 0.001; 
		BC.VonNeumannBC(Boundary::LEFT, *fluid_, 0.001, 0.0);
		BC.VonNeumannBC(Boundary::RIGHT, *fluid_, 0.0, 0.0);

		if (medium_->IsImmersedBodies())
			BC.AdditionalBounceBackBCs();

		//BC.AdditionalBounceBackBCs();

		BC.RecordValuesForAllBC(BCType::PERIODIC, BCType::PERIODIC, BCType::VON_NEUMAN, BCType::VON_NEUMAN);
		if (medium_->IsImmersedBodies())
			BC.RecordAdditionalBCs();
		//BC.RecordAdditionalBCs();

		Recalculate();

		feqCalculate();

		for (auto& i : im_bodies_)
		{
			i->SpreadVelocity(*fluid_);
			i->UpdatePosition();
		}

		std::cout << iter << " Total rho = " << fluid_->rho_.GetSum() << std::endl;

		if (iter % 50 == 0)
		{
			fluid_->WriteFluidToVTK("Data\\ib_lbm_data\\fluid_vtk", iter);

			for (int i = 0; i < im_bodies_.size(); ++i)
			{
				im_bodies_.at(i)->WriteBodyFormToTxt(iter, i);
				im_bodies_.at(i)->WriteBodyFormToVtk("Data\\ib_lbm_data\\body_form_vtk", i, iter);
			}

			fluid_->vx_.WriteFieldToTxt("Data\\ib_lbm_data\\fluid_txt", "vx", iter);
			fluid_->vy_.WriteFieldToTxt("Data\\ib_lbm_data\\fluid_txt", "vy", iter);
			fluid_->rho_.WriteFieldToTxt("Data\\ib_lbm_data\\fluid_txt", "rho", iter);
		}

	}
}

void IBSolver::CreateDataFolder(std::string folder_name) const
{
	// Get path to current directory
	char buffer[MAX_PATH];
	GetModuleFileName(NULL, buffer, MAX_PATH);

	std::string::size_type pos = std::string(buffer).find_last_of("\\/");
	std::string path = std::string(buffer).substr(0, pos);
	path = path.substr(0, path.size() - 6) + "\\" + folder_name;

	char *cstr = new char[path.length() + 1];
	strcpy(cstr, path.c_str());

	// Create folder if not exist yet
	if (GetFileAttributes(cstr) == INVALID_FILE_ATTRIBUTES)
		CreateDirectory(cstr, NULL);
}
