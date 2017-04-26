#include<iostream>
#include<ctime>
#include<cstdlib>
#include<omp.h>

#include"math\2d\my_matrix_2d.h"
#include"phys_values\2d\macroscopic_param_2d.h"
#include"phys_values\2d\distribution_func_2d.h"

#include"modeling_area\medium.h"
#include"modeling_area\fluid.h"
#include"solver\srt.h"
#include"solver\bc\bc.h"


#include"math\3d\my_matrix_3d.h"
#include"phys_values\3d\macroscopic_param_3d.h"
#include"phys_values\3d\distribution_func_3d.h"

int main()
{

#pragma region initial

	srand(time(NULL));

	using std::cout;
	using std::endl;


	omp_set_num_threads(1);
	// ---------------------------

	/*int X{ 100 };
	int Y{ 20 };
	Fluid f(Y, X);
	Medium m(Y, X);

	f.Poiseuille_IC(0.01);

	SRTsolver solver(1.0, m, f);
	solver.solve(100);*/



#pragma endregion

	// 3D matrix testing

	int x{ 6 };
	int y{ 5 };
	int z{ 2 };

	Fluid3D f(z, y, x);
	Medium3D m(z,y,x);

	SRT3DSolver srt(1.0, m, f);
	srt.solve(10);

	/*
		- Продумать структуру для BC!
		- Реализовать ГУ Фон-Неймана для всех границ (TOP, BOTTOM, RIGHT)
		- Не распечатывает файл, т.к. нет дирректории DATA
	*/

	return 0;
}
