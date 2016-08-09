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

	int X{ 15 };
	int Y{ 10 };
	Fluid f(Y, X);
	Medium m(Y, X);

	f.Poiseuille_IC(0.01);

	SRTsolver solver(1.0, m, f);
	solver.solve(50);

	/*
		- Продумать структуру для BC
	*/

	return 0;
}
