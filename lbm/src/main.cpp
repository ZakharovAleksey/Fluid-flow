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

	int x{ 5 };
	int y{ 4 };
	int z{ 3 };


	Medium3D m(z,y,x);
	std::cout << m;


	/*
		- ��������� ��������� ��� BC!
		- ����������� �� ���-������� ��� ���� ������ (TOP, BOTTOM, RIGHT)
		- �� ������������� ����, �.�. ��� ����������� DATA
	*/

	return 0;
}
