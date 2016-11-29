#include<iostream>
#include<ctime>
#include<cstdlib>
#include<omp.h>

#include"math\my_matrix_2d.h"
#include"math\my_matrix_3d.h"

#include"phys_values\macroscopic_param.h"
#include"phys_values\distribution_func.h"
#include"modeling_area\medium.h"
#include"modeling_area\fluid.h"
#include"solver\srt.h"
#include"solver\bc\bc.h"


int main()
{
	srand(time(NULL));

	using std::cout;
	using std::endl;


	omp_set_num_threads(1);

	Matrix3D<int> testMatrix(2, 4, 3);
	testMatrix.Set(1, 3, 2, 10);
	testMatrix(1, 3, 1) = 100;
	std::cout << testMatrix.GetSum() << std::endl;
	


	cout << testMatrix;

	/*
	Matrix2D<int> m(2, 4);
	m(1, 2) = 100;
	cout << m;
	*/
	// 3D tests



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

	/* LBM 2D 
	int X{ 100 };
	int Y{ 20 };
	Fluid f(Y, X);
	Medium m(Y, X);

	f.Poiseuille_IC(0.01);

	SRTsolver solver(1.0, m, f);
	solver.solve(50);
	*/
	

	/*
		- Продумать структуру для BC!
		- Реализовать ГУ Фон-Неймана для всех границ (TOP, BOTTOM, RIGHT)
		- Не распечатывает файл, т.к. нет дирректории DATA
	*/

	return 0;
}
