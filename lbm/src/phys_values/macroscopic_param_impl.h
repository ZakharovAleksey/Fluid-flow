#pragma once

#include"macroscopic_param.h"

template<typename T>
inline MacroscopicParam<T>::MacroscopicParam() : Matrix2D() {}

template<typename T>
MacroscopicParam<T>::MacroscopicParam(unsigned rows, unsigned colls) : Matrix2D(rows, colls) {}
