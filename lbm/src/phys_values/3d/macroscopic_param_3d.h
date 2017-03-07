#pragma once

#ifndef MACROSCOPIC_PARAM_3D_H
#define MACROSCOPIC_PARAM_3D_H

#include"../../math/3d/my_matrix_3d.h"

template<typename T>
class MacroscopicParam3D : public Matrix3D<T>
{
public:
	MacroscopicParam3D();
	MacroscopicParam3D(int depth, int rows, int colls);
	~MacroscopicParam3D();

};

#include"macroscopic_param_3d_impl.h"

#endif // !MACROSCOPIC_PARAM_3D_H
