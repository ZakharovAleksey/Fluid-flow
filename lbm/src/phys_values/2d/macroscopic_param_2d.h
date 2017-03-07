#pragma once

#include"../../math/2d/my_matrix_2d.h"

/*!
	Macroscopic physical values implementation class.

	Implemented on base of Matrix<T> class, i.e. represent Velocity field, or Density field!
*/
template<typename T>
class MacroscopicParam : public Matrix2D<T>
{
public:
	MacroscopicParam();
	MacroscopicParam(unsigned rows, unsigned colls);
	virtual ~MacroscopicParam() {}

private:

};

#include"macroscopic_param_2d_impl.h"