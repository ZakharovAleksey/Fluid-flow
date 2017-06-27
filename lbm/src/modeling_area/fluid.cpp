#include"fluid.h"

#pragma region 2d

Fluid::Fluid(int rows, int colls, const Medium & medium) : rows_(rows), colls_(colls)
{
	rho_.Resize(rows_, colls_);
	FillInitialDensity(medium);

	vx_.Resize(rows_, colls_);
	vy_.Resize(rows_, colls_);

	f_.resize(rows_, colls_);
	feq_.resize(rows_, colls_);
}

Fluid::~Fluid() {}

void Fluid::Poiseuille_IC(double const dvx)
{
	for (int y = 1; y < rows_ - 1; ++y)
		vx_(y, 1) += dvx;
}

void Fluid::WriteFluidToVTK(const std::string path, const int time)
{
	std::string fileName = path + "/fluid_t" + std::to_string(time) + ".vtk";

	std::ofstream output_file;
	output_file.open(fileName);

	if (output_file.is_open())
	{
		// Write VTK header
		// We do not write boundaries, only in [1, size-2] range, because they are equal to zero. So we decrease area size.
		int Nx = GetColumnsNumber() - 2;
		int Ny = GetRowsNumber() - 2;

		output_file << "# vtk DataFile Version 3.0\n";
		output_file << "fluid_state\n";
		output_file << "ASCII\n";
		output_file << "DATASET RECTILINEAR_GRID\n";
		output_file << "DIMENSIONS " << Nx << " " << Ny << " 1" << "\n";

		output_file << "X_COORDINATES " << Nx << " float\n";
		for (int X = 1; X <= Nx; ++X)
			output_file << X + 0.5 << " ";

		output_file << "\n";
		output_file << "Y_COORDINATES " << Ny << " float\n";

		for (int Y = 1; Y <= Ny; ++Y)
			output_file << Y - 0.5 << " ";

		output_file << "\n";
		output_file << "Z_COORDINATES " << 1 << " float\n";
		output_file << 0 << "\n";

		output_file << "POINT_DATA " << Nx * Ny << "\n";
		// Write density
		output_file << "SCALARS Density float 1\n";
		output_file << "LOOKUP_TABLE default\n";
		// On boundaries all macroscopic parameters are equal to zero (this made for good streaming and BC implementation), so 
		// we need to write values in range [1, size - 2]
		for (int Y = Ny; Y >= 1; --Y)
			for (int X = 1; X <= Nx; ++X)
				output_file << rho_(Y, X) << " ";

		output_file << "\nVECTORS Velocity float\n";
		// Write velocity
		for (int Y = Ny; Y >= 1; --Y)
			for (int X = 1; X <= Nx; ++X)
				output_file << vx_(Y, X) << " " << vy_(Y, X) << " 0\n";
	}
	else
	{
		std::cout << "Could not open file for fluid vtk writing!\n";
	}

	output_file.close();
}

void Fluid::FillInitialDensity(const Medium & medium)
{
	// Fills in range [1, size - 1] because of boundaries
	for(int y = 1; y < rows_ - 1; ++y)
		for (int x = 1; x < colls_ - 1; ++x)
		{
			rho_(y, x) = (medium.Get(y, x) != NodeType::BODY_IN_FLUID) ? 1.0 : 0.0;
		}
}

std::pair<unsigned, unsigned> Fluid::size() const
{
	return std::make_pair(rows_, colls_);
}

const int Fluid::GetRowsNumber() const
{
	return rows_;
}

const int Fluid::GetColumnsNumber() const
{
	return colls_;
}


#pragma endregion


#pragma region 3d

Fluid3D::Fluid3D(int depth, int rows, int colls) : depth_(depth), rows_(rows), colls_(colls)
{
	rho_ = std::make_unique<MacroscopicParam3D<double>>(depth_, rows_, colls_);
	rho_->FillWithoutBoundary(1.0);

	vx_ = std::make_unique<MacroscopicParam3D<double>>(depth_, rows_, colls_);
	vy_ = std::make_unique<MacroscopicParam3D<double>>(depth_, rows_, colls_);
	vz_ = std::make_unique<MacroscopicParam3D<double>>(depth_, rows_, colls_);

	f_ = std::make_unique<DistributionFunction3D<double>>(depth_, rows_, colls_);
	feq_ = std::make_unique<DistributionFunction3D<double>>(depth_, rows_, colls_);
}

int Fluid3D::GetDepthNumber() const
{
	return depth_;
}

int Fluid3D::GetRowsNumber() const
{
	return rows_;
}

int Fluid3D::GetColumnsNumber() const
{
	return colls_;
}

void Fluid3D::PoiseuilleIC(double const dvx)
{
	for(int z = 1; z < depth_ - 1; ++ z)
		for (int y = 1; y < rows_ - 1; ++y)
			vy_->operator()(z, y, 1) += dvx;
}

void Fluid3D::SetDistributionFuncValue(const int q, double const value)
{
	assert(q < kQ3d);
	(*f_)[q].FillWith(value);
}

Matrix2D<double> Fluid3D::GetDistributionFuncLayer(const int z, const int q)
{
	Matrix2D<double> res(rows_, colls_);

	for (int y = 0; y < rows_; ++y)
		for (int x = 0; x < colls_; ++x)
			res(y, x) = f_->operator[](q)(z, y, x);
		
	return res;
}

void Fluid3D::SetDistributionFuncLayerValue(const int z, const int q, const int value)
{
	for (int y = 0; y < rows_; ++y)
		for (int x = 0; x < colls_; ++x)
			(*f_)[q](z, y, x) = value;
}

void Fluid3D::RecalculateRho()
{
	rho_->FillWith(0.0);
	for (int q = 0; q < kQ3d; ++q)
	{
		*rho_ += (*f_)[q];
	}
}

void Fluid3D::RecalculateV()
{
	RecalculateVelocityComponent(vx_, ex);
	RecalculateVelocityComponent(vy_, ey);
	RecalculateVelocityComponent(vz_, ez);
}

long double Fluid3D::TotalRho()
{
	return rho_->GetSum();
}

void Fluid3D::WriteFluidToVTK(const std::string path, int time)
{
	// Create filename
	std::string fileName = path + "/fluid_t" + std::to_string(time) + ".vtk";

	std::ofstream output_file;
	output_file.open(fileName);

	if (output_file.is_open())
	{
		// Write VTK header

		// We do not write boundaries, only in [1, size-2] range, because they are equal to zero. So we decrease area size.
		const int Nx = GetRowsNumber() - 2;
		const int Ny = GetColumnsNumber() - 2;
		const int Nz = GetDepthNumber() - 2;

		output_file << "# vtk DataFile Version 3.0\n";
		output_file << "3d_srt_fluid_state\n";
		output_file << "ASCII\n";
		output_file << "DATASET RECTILINEAR_GRID\n";
		output_file << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";

		output_file << "X_COORDINATES " << Nx << " float\n";
		for (int X = 0; X < Nx; ++X)
			output_file << X + 0.5 << " ";
		output_file << "\n";

		output_file << "Y_COORDINATES " << Ny << " float\n";
		for (int Y = 0; Y < Ny; ++Y)
			output_file << Y + 0.5 << " ";
		output_file << "\n";

		output_file << "Z_COORDINATES " << Nz << " float\n";
		for (int Z = 0; Z < Nz; ++Z)
			output_file << Z + 0.5 << " ";
		output_file << "\n";


		output_file << "POINT_DATA " << Nx * Ny * Nz << "\n";
		// Write Density(scalar)
		output_file << "SCALARS Density float 1\n";
		output_file << "LOOKUP_TABLE default\n";
		// On boundaries all macroscopic parameters are equal to zero (this made for good streaming and BC implementation), so 
		// we need to write values in range [1, size - 2]
		for (int Z = Nz; Z > 0; --Z)
			for (int Y = Ny; Y > 0; --Y)
				for (int X = 1; X <= Nx; ++X)
				{
					output_file << rho_->operator()(Z, Y, X) << " ";
				}

		// Write Velocity(vector V = (Vx, Vy, Vz))
		output_file << "\nVECTORS Velocity float\n";
		for (int Z = Nz; Z > 0; --Z)
			for (int Y = Ny; Y > 0; --Y)
				for (int X = 1; X <= Nx; ++X)
					output_file << vx_->operator()(Z, Y, X) << " " << vy_->operator()(Z, Y, X) << " " << vz_->operator()(Z, Y, X) << "\n";
	}
	else
	{
		std::cout << "Could not open file for fluid vtk writing!\n";
	}

	output_file.close();
}
void Fluid3D::RecalculateVelocityComponent(const MacroscopicParamPtr & v_ptr, const int e[])
{
	v_ptr->FillWith(0.0);

	for (int q = 0; q < kQ3d; ++q)
		*v_ptr += (*f_)[q] * (double) e[q];

	v_ptr->TimesDivide(*rho_);

}


#pragma endregion

