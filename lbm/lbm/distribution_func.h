#pragma once

#include<array>

#include"my_matrix.h"


/*!
	Возможное число направлений, в которых могут перемещаться
	псевдочастицы. Следует из размерности задачи D2Q3 D2Q9.
*/
unsigned const kQ{ 9 };

/*!
	Функция распределения плотности вероятности f.
		- массив размера kQ в каждой ячейке которого хранятся значения функции распределения для q = 0...kQ
		  направления.
*/
template<typename T>
class DistributionFunction
{
public:
	DistributionFunction();
	DistributionFunction(unsigned rows, unsigned colls);
	~DistributionFunction();
	
	DistributionFunction(DistributionFunction<T> const & other);
	
	void swap(DistributionFunction & dist_func);

	DistributionFunction<T> & operator=(DistributionFunction<T> const & other) {
		assert(rows_ == other.rows_ && colls_ == other.colls_);
		if (this != &other) {
			DistributionFunction<T> temp(other);
			temp.swap(*this);
		}
		return *this;
	}

	DistributionFunction<T> & operator+=(DistributionFunction<T> const & other) {
	#pragma omp parallel for
		for (int q = 0; q < kQ; ++q)
			dfunc_body_.at(q) += other.dfunc_body_.at(q);

		return *this;
	}

	DistributionFunction<T> operator+(DistributionFunction<T> const & other) {
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		DistributionFunction<T> result(*this);
	#pragma omp parallel for
		for (int q = 0; q < kQ; ++q)
			result.dfunc_body_.at(q) += other.dfunc_body_.at(q);

		return result;
	}


	DistributionFunction<T> & operator-=(DistributionFunction<T> const & other) {
	#pragma omp parallel for
		for (int q = 0; q < kQ; ++q)
			dfunc_body_.at(q) -= other.dfunc_body_.at(q);

		return *this;
	}

	DistributionFunction<T> operator-(DistributionFunction<T> const & other) {
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		DistributionFunction<T> result(*this);
	#pragma omp parallel for
		for (int q = 0; q < kQ; ++q)
			result.dfunc_body_.at(q) -= other.dfunc_body_.at(q);

		return result;
	}


	DistributionFunction<T> & operator/=(T const & other) {
	#pragma omp parallel for
		for (int q = 0; q < kQ; ++q)
			dfunc_body_.at(q) /= other;

		return *this;
	}

	DistributionFunction<T> operator/(T const & other) {
		
		DistributionFunction<T> result(*this);
	#pragma omp parallel for
		for (int q = 0; q < kQ; ++q)
			result.dfunc_body_.at(q) /= other;

		return result;
	}

	//! Возвращвет q-ую компоненту функции распределения
	Matrix<T> & operator[](unsigned q);

	//! Заполняетя все элементы массивов функций распределения кроме элементов на границе величиной value
	void fill(T const value);

	//! Изменяет размер каждой из компонент функции распределения на rows и colls соответственно
	void resize(unsigned rows, unsigned colls);


	/*! 
		В q-yю строку возвращаемой матрицы записывает величины q-ой компоненты функции распределения 
		находящейся на верхней границе.
	*/
	std::array<Matrix<T>, kQ> get_values_on_upper_boundary() const;

	/*!
	В q-yю строку возвращаемой матрицы записывает величины q-ой компоненты функции распределения
	находящейся на нижней границе.
	*/
	std::array<Matrix<T>, kQ> get_values_on_bottom_boundary() const;

	/*!
	В q-yю строку возвращаемой матрицы записывает величины q-ой компоненты функции распределения
	находящейся на левой границе.
	*/
	std::array<Matrix<T>, kQ> get_values_on_left_boundary() const;

	/*!
	В q-yю строку возвращаемой матрицы записывает величины q-ой компоненты функции распределения
	находящейся на правой границе.
	*/
	std::array<Matrix<T>, kQ> get_values_on_right_boundary() const;


	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, DistributionFunction<T1> const & distr_func);

private:
	unsigned rows_;
	unsigned colls_;
	
	/*!
		Хранит значения функции распределения плотности вероятности вдоль q = 1...kQ направления движения
		псевдочастицы.

		Реализован как массив матриц фиксированного рамера kQ.
	*/
	std::array<Matrix<T>, kQ> dfunc_body_;
};


#include"distribution_func_impl.h"