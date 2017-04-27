#include"fluid.h"

#pragma region 2d

Fluid::Fluid(unsigned rows, unsigned colls) : rows_(rows), colls_(colls) {

	rho_.Resize(rows_, colls_);
	rho_.FillWithoughtBoundary(1.0);
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

std::pair<unsigned, unsigned> Fluid::size() const
{
	return std::make_pair(rows_, colls_);
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

void Fluid3D::Poiseuille_IC(double const dvx)
{
	for(int z = 1; z < depth_ - 1; ++ z)
		for (int y = 1; y < rows_ - 1; ++y)
			vx_->operator()(z, y, 1) += dvx;
}

void Fluid3D::SetDistributionFuncValue(const int q, double const value)
{
	assert(q < kQ3d);
	f_->operator[](q).FillWith(value);
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
			f_->operator[](q)(z, y, x) = value;
}

void Fluid3D::CalculateDencity()
{
	for(int z = 0; z < depth_; ++z)
		for(int y = 0; y < rows_; ++y)
			for (int x = 0; x < colls_; ++x)
			{
				rho_->operator()(z, y, x) = 0.0;
				for (int q = 0; q < kQ3d; ++q)
					rho_->operator()(z, y, x) += f_->operator[](q)(z, y, x);
			}
}

void Fluid3D::CalculateVelocity(const MacroscopicParam3D<double>& rho, const int[] e)
{

	for (int q = 0; q < kQ; ++q)
		result += dfunc_body_.at(q) * mas[q];

	// Переписать сразу в return, тк лишнее копирование НО выпадает error
	result.TimesDivide(density);
	return result;
}


#pragma endregion

