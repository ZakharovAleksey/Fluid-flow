#include<iostream>
#include<ctime>
#include<cstdlib>

#include<omp.h>

#include"my_matrix.h"

#include"macroscopic_param.h"

#include"distribution_func.h"

#include"medium.h"

#include"fluid.h"

#include"srt.h"

#include"bc.h"

int main()
{

	srand(time(NULL));
	using std::cout;
	using std::endl;

	omp_set_num_threads(1);

	// Тестирование функционала матриц ------ 

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

	// Distribution function test ------- 

	/*DistributionFunction<double> dfunc(5, 10);
	cout << dfunc;

	std::vector<double> a(dfunc.get_top_boundary_val(1));
	for (auto i : a)
		cout << i << ' ';
	cout << endl;

	std::vector<double> b(dfunc.get_bottom_boundary_val(1));
	for (auto i : b)
		cout << i << ' ';
	cout << endl;

	std::vector<double> c(dfunc.get_right_boundary_val(1));
	for (auto i : c)
		cout << i << ' ';
	cout << endl;

	std::vector<double> d(dfunc.get_left_boundary_val(1));
	for (auto i : d)
		cout << i << ' ';
	cout << endl;*/

	// ----------------------------

	int X{ 100 };
	int Y{ 20 };
	Fluid f(Y, X);
	Medium m(Y, X);

	f.Poiseuille_IC(0.01);

	SRTsolver solver(1.0, m, f);
	solver.solve(50);

	/*
		- Продумать структуру для BC
		- Тестировать на задаче пуазейля! Если работет нормально то тогда устанавливать библиотеки к Python
		
		- Возможно погрешности в плотности из-за того, что у нас в Фон-Неймане скорость пересчитывается а плотность нет
	*/

	return 0;
}
