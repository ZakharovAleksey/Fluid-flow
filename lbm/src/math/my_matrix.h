#pragma once

#include<iostream>
#include<ostream>
#include <cassert> // setw()
#include<vector>
#include <functional>
#include<string>
#include<omp.h>

#define TEST_INCLUDE 

#ifndef TEST_INCLUDE 
 	#include<cstdlib>
 	#include<ctime>
#endif // !TEST_INCLUDE


/*!
	Класс матрицы на базе которой будут реализованны физические величины, такие как плотность, скорость, 
	функция распределения плотности вероятности.

	Перегрузка производилась ТОЛЬКО тех операторов, котрые применяются в методе.

*/
template<typename T>
class Matrix
{
public:
	Matrix();
	Matrix(unsigned rows, unsigned colls);
	virtual ~Matrix();

	Matrix(Matrix<T> const & other);	// Проверить все ли хорошо при вызове функции отсюда

	void swap(Matrix<T> & other);

#pragma region operator+

	Matrix<T> & operator=(Matrix<T> const & other) 
	{
		assert(rows_ == other.rows_ && colls_ == other.colls_);
		if (this != & other) {
			Matrix<T> temp(other);
			temp.swap(*this);
		}

		return *this;
	}

	Matrix<T> & operator+=(Matrix<T> const & other) {
		assert(rows_ == other.rows_ && colls_ == other.colls_);
	#pragma omp parallel for
			for (int i = 0; i < body_.size(); ++i)
				body_.at(i) += other.body_.at(i);

		return *this;
	}

	Matrix<T> & operator+=(T const other) {
		std::for_each(body_.begin(), body_.end(),
			[&](T & value) {value += other; });

		return *this;
	}

	Matrix<T> operator+(Matrix<T> const & other) const {
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		Matrix<T> result(*this);
	#pragma omp parallel for
		for (int i = 0; i < result.body_.size(); ++i)
			result.body_.at(i) += other.body_.at(i);

		return result;
	}

	Matrix<T> operator+(T const other) const {
		Matrix<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value += other; });

		return result;
	}

	template<typename T1>
	friend Matrix<T1> operator+(T1 const left, Matrix<T1> const & right);

#pragma endregion

#pragma region operator-

	Matrix<T> & operator-=(Matrix<T> const & other) {
		assert(rows_ == other.rows_ && colls_ == other.colls_);
	#pragma omp parallel for
			for (int i = 0; i < body_.size(); ++i)
				body_.at(i) -= other.body_.at(i);

		return *this;
	}

	Matrix<T> & operator-=(T const other) {
		std::for_each(body_.begin(), body_.end(),
			[&](T & value) {value -= other; });

		return *this;
	}

	Matrix<T> operator-(Matrix<T> const & other) const {
		assert(rows_ == other.rows_ && colls_ == other.colls_);

		Matrix<T> result(*this);
	#pragma omp parallel for
		for (int i = 0; i < result.body_.size(); ++i)
			result.body_.at(i) -= other.body_.at(i);

		return result;
	}

	Matrix<T> operator-(T const other) const {
		Matrix<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value -= other; });

		return result;
	}

#pragma endregion

#pragma region operator*

	Matrix<T> & operator*=(T other) {
	#pragma omp parallel for
		for (int i = 0; i < body_.size(); ++i)
			body_.at(i) *= other;

		return *this;
	}

	Matrix<T> operator*(T other) const {
		Matrix<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value *= other; });

		return result;
	}

	template<typename T1>
	friend Matrix<T1> operator*(T1 const left, Matrix<T1> const & right);

	/*!
	Скалярное произведение двух матриц то есть соответсвующий элемент первой матрицы
	умножается на соответствующий элемент второй матрицы.
	*/
	Matrix<T> scalar_multiplication(Matrix<T> const & other);

#pragma endregion

#pragma region operator/

	Matrix<T> & operator/=(T other) {
	#pragma omp parallel for
		for (int i = 0; i < body_.size(); ++i)
			body_.at(i) /= other;

		return *this;
	}

	Matrix<T> operator/(T other) const {
		Matrix<T> result(*this);

		std::for_each(result.body_.begin(), result.body_.end(),
			[&](T & value) {value /= other; });

		return result;
	}

	//! Почленное деление на матрицу в аргументе
	Matrix<T> times_divide(Matrix<T> const & other);

#pragma endregion

	T & operator()(unsigned y, unsigned x) {
		return body_.at(y * colls_ + x);
	}

	T const operator()(unsigned y, unsigned x) const {
		return body_.at(y * colls_ + x);
	}

#pragma region functions

	//! Возвращает пару в которой first = rows_, second = colls_
	std::pair<unsigned int, unsigned int> size() const;

	//! Возвращает сумму элементов матрицы
	long double get_sum() const;

	//! Возвращает std::vector<T>, в которой хранятся все элементы строки с номером y
	std::vector<T> get_row(unsigned const y) const;

	//! Задает элементы матрицы строки y равными элементам std::vector<T> row
	void set_row(unsigned const y, std::vector<T> const & row);

	//! Возвращает std::vector<T>, в которой хранятся все элементы диапазона [1 : rows_ - 2] столбца с номером x
	std::vector<T> get_coll(unsigned const x) const;

	//! Задает элементы диапазона [1 : colls - 1] матрицы столбца x равными элементам std::vector<T> coll
	void set_coll(unsigned const x, std::vector<T> const & coll);

	//! Заполняет матрицу целиком значением value
	void fill_with(T const value);

	//! Заполняет матрицу без граничных элементов значениями value
	void fill_withought_boundary(T const value);

	//! Заполняет столбец номер coll_id матрицы значениями value
	void fill_coll_with(int const coll_id, T const value);

	//! Заполняет строку номер row_id матрицы значениями value
	void fill_row_with(int const row_id, T const value);

	//! Изменяет размеры матрицы на rows и colls соответственно, заполняя при этом нулями
	void resize(unsigned rows, unsigned colls);

#pragma endregion

#pragma ostream_functions

	//! Выводит объект типа матрица в файл
	void to_file(std::string value_name, int const time);

	//! Печатает столбец номер coll_id в файл
	void coll_to_file(std::string value_name, int const coll_id, int const time);

	//! Печатает столбец номер coll_id в файл
	void row_to_file(std::string value_name, int const row_id, int const time);

#pragma endregion

	template<typename T1>
	friend std::ostream & operator<<(std::ostream & os, Matrix<T1> const & matrix);

protected:
	unsigned rows_;
	unsigned colls_;

	std::vector<T> body_;
};

#pragma region additional_functions
	
#pragma endregion


#include"my_matrix_impl.h"
