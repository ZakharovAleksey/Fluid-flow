#include<iostream>
#include<ctime>
#include<cstdlib>
#include<omp.h>

#include"math\2d\my_matrix_2d.h"
#include"phys_values\macroscopic_param.h"
#include"phys_values\distribution_func.h"
#include"modeling_area\medium.h"
#include"modeling_area\fluid.h"
#include"solver\srt.h"
#include"solver\bc\bc.h"


#include"math\3d\my_matrix_3d.h"

int main()
{
	srand(time(NULL));

	using std::cout;
	using std::endl;


	omp_set_num_threads(5);

	// Тестирование функционала матриц ------ Реализовать через тесты

	/*Matrix<double> m(5, 4);
	cout << m;
	cout << "M rows = " << m.size().first << " M colls = " << m.size().second << endl;
	std::vector<double> line = m.get_row(2);
	for (auto i : line)
		cout << i;
	std::cout << std::endl;
	std::vector<double> row = m.get_coll(2);
	for (auto i : row)
		cout << i;
	std::cout << std::endl;
	m.set_coll(3, row);
	m.set_row(0, line);
	cout << m;*/

	// ---------------------------

	/*int X{ 100 };
	int Y{ 20 };
	Fluid f(Y, X);
	Medium m(Y, X);

	f.Poiseuille_IC(0.01);

	SRTsolver solver(1.0, m, f);
	solver.solve(10);*/


	// 3D matrix testing

	int x{ 3 };
	int y{ 6 };
	int z{ 3 };

	Matrix3D<int> first(z, y, x);
	Matrix3D<int> second(x, y, z);

	std::cout << first;
	std::vector<int> row = first.GetRow(2);
	std::vector<int> column = first.GetColumn(2);

	for (auto & i : column)
	{
		std::cout << i << " ";
		i += 100;
	}
	std::cout << std::endl;

	first.SetColumn(1, column);

	std::cout << first;



	/*for (auto i : row)
		std::cout << i << " ";
	std::cout << std::endl;

	for (auto & i : row)
		i += 100;

	first.SetRow(3, row);

	std::cout << first;*/

	/*
		- Продумать структуру для BC!
		- Реализовать ГУ Фон-Неймана для всех границ (TOP, BOTTOM, RIGHT)
		- Не распечатывает файл, т.к. нет дирректории DATA
	*/

	return 0;
}
