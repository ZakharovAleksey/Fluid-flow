#pragma once

#ifndef MACROSCOPIC_PARAM_3D_IMPL_H
#define MACROSCOPIC_PARAM_3D_IMPL_H

#include"macroscopic_param_3d.h"

template<typename T>
MacroscopicParam3D<T>::MacroscopicParam3D() : Matrix3D<T>() {}

template<typename T>
MacroscopicParam3D<T>::MacroscopicParam3D(int depth, int rows, int colls) : Matrix3D<T>(depth, rows, colls) {}

template<typename T>
MacroscopicParam3D<T>::~MacroscopicParam3D() {}

#endif // !MACROSCOPIC_PARAM_3D_IMPL_H
