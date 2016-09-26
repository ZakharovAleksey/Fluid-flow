#pragma once

#include"..\math\my_matrix.h"

/*!
	Macroscopic physical values implementation class.

	Implemented on base of Matrix<T> class, i.e. represent Velocity field, or Density field!
*/
template<typename T>
class MacroscopicParam : public Matrix<T>
{
public:
	MacroscopicParam();
	MacroscopicParam(unsigned rows, unsigned colls);
	virtual ~MacroscopicParam() {}

private:

};

#include"macroscopic_param_impl.h"