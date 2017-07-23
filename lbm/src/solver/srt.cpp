#include"srt.h"


#pragma region 2d


#pragma region srt

SRTsolver::SRTsolver(double const tau, Medium & medium, Fluid & fluid, BCs* bc) : tau_(tau), medium_(&medium), fluid_(&fluid), bc_(bc)
{
	assert(medium_->size().first == fluid_->size().first);
	assert(medium_->size().second == fluid_->size().second);

	CreateDataFolder("Data");
	CreateDataFolder("Data\\srt_lbm_data");
	CreateDataFolder("Data\\srt_lbm_data\\2d");
	CreateDataFolder("Data\\srt_lbm_data\\2d\\fluid_txt");
	CreateDataFolder("Data\\srt_lbm_data\\2d\\fluid_vtk");
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

void SRTsolver::Streaming()
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

	// Очищаем значения попавшие на границу, так как они уже сохранены в BCs
	fluid_->f_.fillBoundaries(0.0);
}

void SRTsolver::Collision()
{
	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] += (fluid_->feq_[q] - fluid_->f_[q]) / tau_;
}

void SRTsolver::Solve(int iter_numb)
{
	feqCalculate();

	for (int q = 0; q < kQ; ++q)
		fluid_->f_[q] = fluid_->feq_[q];

	std::cout << " Total rho = " << fluid_->rho_.GetSum() << std::endl;

	for (int iter = 0; iter < iter_numb; ++iter) 
	{
		Collision();

		bc_->PrepareBCValues(*medium_);
		Streaming();
		bc_->PerformBC(fluid_, *medium_);
		bc_->RecordBCValues(*medium_);

		Recalculate();
		
		feqCalculate();

		std::cout << iter << " Total rho = " << fluid_->rho_.GetSum() << std::endl;

		if (iter % 100 == 0)
		{
			fluid_->vx_.WriteFieldToTxt("Data\\srt_lbm_data\\2d\\fluid_txt", "vx", iter);
			fluid_->vy_.WriteFieldToTxt("Data\\srt_lbm_data\\2d\\fluid_txt", "vy", iter);
			fluid_->WriteFluidToVTK("Data\\srt_lbm_data\\2d\\fluid_vtk", iter);
		}

	}
}

void SRTsolver::Recalculate()
{
	fluid_->rho_ = fluid_->f_.calculateDensity();
	fluid_->vx_ = fluid_->f_.calculateVelocity(kEx, fluid_->rho_);
	fluid_->vy_ = fluid_->f_.calculateVelocity(kEy, fluid_->rho_);
}

void SRTsolver::CreateDataFolder(std::string folder_name) const
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

#pragma endregion




#pragma endregion

#pragma region 3d


SRT3DSolver::SRT3DSolver(double tau, Medium3D & medium, Fluid3D & fluid) : tau_(tau), medium_(& medium), fluid_(& fluid)
{
	assert(medium_->GetDepthNumber() == fluid_->GetDepthNumber());
	assert(medium_->GetRowsNumber() == fluid_->GetRowsNumber());
	assert(medium_->GetColumnsNumber() == fluid_->GetColumnsNumber());

	CreateDataFolder("Data");
	CreateDataFolder("Data\\srt_lbm_data");
	CreateDataFolder("Data\\srt_lbm_data\\3d");
	CreateDataFolder("Data\\srt_lbm_data\\3d\\fluid_txt");
	CreateDataFolder("Data\\srt_lbm_data\\3d\\fluid_vtk");
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

void SRT3DSolver::Streaming()
{
	const int depth = medium_->GetDepthNumber();
	const int rows = medium_->GetRowsNumber();
	const int colls = medium_->GetColumnsNumber();

	SubStreamingMiddle(depth, rows, colls);
	SubStreamingTop(depth, rows, colls);
	SubStreamingBottom(depth, rows, colls);
}

void SRT3DSolver::Collision()
{
	for (int q = 0; q < kQ3d; ++q)
		(*fluid_->f_)[q] += ((*fluid_->feq_)[q] - (*fluid_->f_)[q]) / tau_;
}

void SRT3DSolver::Solve(int iter_numb)
{
	//fluid_->vz_->SetTBLayer(1, std::vector<double>(fluid_->GetColumnsNumber() * fluid_->GetRowsNumber(), 0.01));

	feqCalculate();
	for (int q = 0; q < kQ; ++q)
		(*fluid_->f_)[q] = (*fluid_->feq_)[q];

	BCs3D bc(fluid_->GetRowsNumber(), fluid_->GetColumnsNumber(), *fluid_->f_);

	for (int iter = 0; iter < iter_numb; ++iter)
	{
		std::cout << iter << " : ";
		Collision();
		bc.PrepareValuesForAllBC(BCType::VON_NEUMAN, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK);

		Streaming();

		//bc.BounceBackBC(Boundary::TOP);
		bc.VonNeumannBC(Boundary::TOP, 0.0, 0.0, 0.001);
		bc.BounceBackBC(Boundary::BOTTOM);
		bc.BounceBackBC(Boundary::LEFT);
		bc.BounceBackBC(Boundary::RIGHT);
		bc.BounceBackBC(Boundary::CLOSE_IN);
		bc.BounceBackBC(Boundary::FAAR);

		bc.RecordValuesForAllBC(BCType::VON_NEUMAN, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK, BCType::BOUNCE_BACK);

		Recalculate();
		//fluid_->vz_->SetTBLayer(1, std::vector<double>(fluid_->GetColumnsNumber() * fluid_->GetRowsNumber(), 0.01));
		

		feqCalculate();

		if (iter % 10 == 0)
			fluid_->WriteFluidToVTK("Data\\srt_lbm_data\\3d\\fluid_vtk", iter);
		if(iter % 50 == 0)
			GetNFProfile(*fluid_->vz_, 5, iter);
	}
	
}

void SRT3DSolver::GetProfile(const int chan_numb, const int iter_numb)
{
	std::vector<double> res = fluid_->vy_->GetNFLayer(chan_numb);
	int colls = fluid_->GetColumnsNumber();
	int rows = fluid_->GetRowsNumber();
	int depth = fluid_->GetDepthNumber();

	

	std::string name = "Data\\srt_lbm_data\\3d\\fluid_txt\\profile_" + std::to_string(iter_numb) + ".txt";
	if (WriteHeatMapInFile(name, res, colls - 2))
	{
		std::cout << "Data writing complete successfully!\n";
	}
	
	
}

void SRT3DSolver::GetNFProfile(const MacroscopicParam3D<double> & physVal, const int chan_numb, const int iter_numb)
{
	std::vector<double> res = physVal.GetNFLayer(chan_numb);

	int colls = fluid_->GetColumnsNumber();
	int rows = fluid_->GetRowsNumber();
	int depth = fluid_->GetDepthNumber();

	std::string name = "Data\\srt_lbm_data\\3d\\fluid_txt\\profile_tb_" + std::to_string(iter_numb) + ".txt";
	if (WriteHeatMapInFile(name, res, colls - 2))
	{
		std::cout << "Data writing complete successfully!\n";
	}


}

bool SRT3DSolver::WriteHeatMapInFile(const std::string & file_name, const std::vector<double>& data, const int lenght)
{
	std::ofstream file;
	file.open(file_name);
	if (file.is_open())
	{
		file.precision(3);

		for (int i = 0; i < data.size(); ++i)
		{
			// If this is laset element in the row
			if (i != 0 && ((i + 1) % lenght == 0))
			{
				// If it is not last element in list we need to and endl
				if (i != data.size() - 1)
					file << data.at(i) << std::endl;
				// If it is las element in list we need no endl
				else
					file << data.at(i);
			}
			else
				file << data.at(i) << " ";
		}

		file.close();
		return true;
	}
	else
	{
		std::cout << "Could not open file: " << file_name << " on writing!\n";
	}
	return false;

	
}

void SRT3DSolver::SubStreamingMiddle(const int depth, const int rows, const int colls)
{
	for (int z = 0; z < depth; ++z)
	{
		const std::vector<int> middle_ids_{ 0,1,2,3,6,7,10,11,18};
		for (auto q : middle_ids_)
		//for (int q = 0; q < 9; ++q)
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
	//for (int z = 1; z < depth; ++z)
	for (int z = depth - 2; z > 0 ; --z)
	{
		const std::vector<int> top_ids_{ 4,8,12,14,17 };
		for(auto q : top_ids_)
		//for (int q = 9; q < 14; ++q)
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
	//for (int z = depth - 1; z > 0; --z)
	for (int z = 1; z < depth - 1; ++z)
	{
		const std::vector<int> bottom_ids_{ 5,9,13,15,16}; 
		for (auto q : bottom_ids_)
		//for (int q = 14; q < 19; ++q)
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

void SRT3DSolver::CreateDataFolder(std::string folder_name) const
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

void SRT3DSolver::Recalculate()
{
	fluid_->RecalculateRho();
	fluid_->rho_->FillTopBottomWalls(0.0);
	fluid_->RecalculateV();

	std::cout << "Total rho: " << fluid_->TotalRho() << std::endl;
}

#pragma endregion
