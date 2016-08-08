#pragma once

#include"macroscopic_param.h"

template<typename T>
inline MacroscopicParam<T>::MacroscopicParam() : Matrix() {}

template<typename T>
MacroscopicParam<T>::MacroscopicParam(unsigned rows, unsigned colls) : Matrix(rows, colls) {}
