#pragma once

#include"macroscopic_param_2d.h"

template<typename T>
inline MacroscopicParam<T>::MacroscopicParam() : Matrix2D() {}

template<typename T>
MacroscopicParam<T>::MacroscopicParam(unsigned rows, unsigned colls) : Matrix2D(rows, colls) {}
