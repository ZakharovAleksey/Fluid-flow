#include"fluid.h"

Fluid::Fluid(unsigned rows, unsigned colls) : rows_(rows), colls_(colls) {

	rho_.resize(rows_, colls_);
	rho_.fill_withought_boundary(1.0);
	vx_.resize(rows_, colls_);
	vy_.resize(rows_, colls_);

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
