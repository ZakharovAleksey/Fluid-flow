#pragma once

#include<array>

#include"..\math\my_matrix.h"
#include"..\phys_values\macroscopic_param.h"

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

#pragma region opeartors

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

#pragma endregion

#pragma region function

	//! Возвращвет q-ую компоненту функции распределения
	Matrix<T> & operator[](unsigned q);

	//! Заполняетя все элементы массивов функций распределения кроме элементов на границе величиной value
	void fill(T const value);

	//! Заполняет границы для каждой из kQ функций распределения значениями value
	void fill_boundaries(T const value);

	//! Изменяет размер каждой из компонент функции распределения на rows и colls соответственно
	void resize(unsigned rows, unsigned colls);

	//! Возвращает пару в которой first = rows_, second = colls_
	std::pair<unsigned int, unsigned int> size() const;

	//! Возвращает значения функции распределения на верхней стенке
	std::vector<T> get_top_boundary_val(int const q) const;	
	//! Возвращает значения функции распределения на нижней стенке
	std::vector<T> get_bottom_boundary_val(int const q) const;
	//! Возвращает значения функции распределения на левой стенке
	std::vector<T> get_left_boundary_val(int const q) const;
	//! Возвращает значения функции распределения на правой стенке
	std::vector<T> get_right_boundary_val(int const q) const;
	
	//! Назначет значения функции распределения на верхней стенке равными значению масива colls
	void set_top_boundary_value(int const q, std::vector<T> const & row);
	//! Назначет значения функции распределения на нижней стенке равными значению масива colls
	void set_bottom_boundary_value(int const q, std::vector<T> const & row);
	//! Назначет значения функции распределения на левой стенке равными значению масива colls
	void set_left_boundary_value(int const q, std::vector<T> const & coll);
	//! Назначет значения функции распределения на правой стенке равными значению масива colls
	void set_right_boundary_value(int const q, std::vector<T> const & coll);

	//! Считает плотность в кажой из ячеек области
	MacroscopicParam<T> get_density() const;
	//! Считает скорость в кажой из ячеек области
	MacroscopicParam<T> get_velocity(const double mas[kQ], MacroscopicParam<T> const & density) const;

#pragma endregion

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