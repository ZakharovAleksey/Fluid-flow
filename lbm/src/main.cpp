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

	int x{ 2 };
	int y{ 2 };
	int z{ 2 };

	Matrix3D<int> first(x, y, z);
	Matrix3D<int> second(x, y, z);

	std::cout << first;
	std::cout << second;
	
	first = 10 + 10* second / 10 - 10;


	std::cout << first;

	/*
		- Продумать структуру для BC!
		- Реализовать ГУ Фон-Неймана для всех границ (TOP, BOTTOM, RIGHT)
		- Не распечатывает файл, т.к. нет дирректории DATA
	*/

	return 0;
}
