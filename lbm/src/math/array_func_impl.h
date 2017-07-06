#pragma once

#include <algorithm>
#include <functional>

/*
	Перегружает элементарные операции для std::vector<>
	 - operator+
	 - operator-
	 - operator*
	 - operator/
*/

template <typename T>
std::vector<T> operator+(std::vector<T> const & left, std::vector<T> const & right)
{
	assert(left.size() == right.size());

	std::vector<T> result(left);

	for (int i = 0; i < result.size(); ++i)
		result.at(i) += right.at(i);

	return result;
}

template <typename T>
std::vector<T> operator+(T const & left, std::vector<T> const & right)
{
	std::vector<T> result(right);

	for (int i = 0; i < result.size(); ++i)
		result.at(i) += left;

	return result;
}

template <typename T>
std::vector<T> operator-(T const & left, std::vector<T> const & right)
{

	std::vector<T> result(right.size(), left);

	for (int i = 0; i < result.size(); ++i)
		result.at(i) -= right.at(i);

	return result;
}

template<typename T>
std::vector<T>& operator+=(Matrix3D<T>& left, const Matrix3D<T>& right)
{
	assert(left.size() == right.size());

	for (int i = 0; i < result.size(); ++i)
		left.at(i) += right.at(i);

	return left;
}


template <typename T>
std::vector<T> operator-(std::vector<T> const & left, std::vector<T> const & right)
{
	assert(left.size() == right.size());

	std::vector<T> result(left);

	for (int i = 0; i < result.size(); ++i)
		result.at(i) -= right.at(i);

	return result;
}

template <typename T>
std::vector<T> operator-(std::vector<T> const & left, T const right)
{
	std::vector<T> result(left);

	for (int i = 0; i < left.size(); ++i)
		result.at(i) -= right;

	return result;
}

template <typename T>
std::vector<T> operator+(std::vector<T> const & left, T const right)
{
	std::vector<T> result(left);

	for (int i = 0; i < left.size(); ++i)
		result.at(i) += right;

	return result;
}

template <typename T>
std::vector<T> operator*(std::vector<T> const & left, std::vector<T> const & right)
{
	std::vector<T> result(right);

	for (int i = 0; i < right.size(); ++i)
		result.at(i) *= left.at(i);
	return result;
}

template <typename T>
std::vector<T> operator*(T const left, std::vector<T> const & right)
{
	std::vector<T> result(right);

	for (int i = 0; i < right.size(); ++i)
		result.at(i) *= left;
	return result;
}

template <typename T>
std::vector<T> operator*(std::vector<T> const & left, T const right)
{
	std::vector<T> result(left);

	for (int i = 0; i < left.size(); ++i)
		result.at(i) *= right;

	return result;
}

template <typename T>
std::vector<T> operator/(std::vector<T> const & left, T const right)
{
	std::vector<T> result(left);

	for (int i = 0; i < left.size(); ++i)
		result.at(i) /= right;
	return result;
}

template <typename T>
std::vector<T> operator/(std::vector<T> const & left, std::vector<T> const & right)
{
	std::vector<T> result(left);

	for (int i = 0; i < left.size(); ++i)
		result.at(i) /= right.at(i);
	return result;
}

//! Выводит значения вектора на экран
template<typename T>
void  display(std::ostream & os, std::vector<T> const & body) {
	os <<"Size = " << body.size() << std::endl;
	std::copy(body.begin(), body.end(), std::ostream_iterator<T>(os, " "));
}
