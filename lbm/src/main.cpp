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



void MatrixTest()
{
	int x{ 6 };
	int y{ 5 };
	int z{ 4 };

	Matrix3D<int> m(z, y, x);

	for (int zz = 0; zz < z; ++zz)
		for (int yy = 0; yy < y; ++yy)
			for (int xx = 0; xx < x; ++xx)
				m(zz, yy, xx) = rand() % 100;

	std::cout << m << std::endl;


	std::vector<int> layer = m.GetTBLayer(z - 2);
	for (auto & i : layer)
	{
		std::cout << i << " ";
		i = 0;
	}
	std::cout << std::endl;

	std::vector<int> layer1 = m.GetLRLayer(x - 2);
	for (auto & i : layer1)
	{
		std::cout << i << " ";
		i = 1;
	}
	std::cout << std::endl;

	std::vector<int> layer2 = m.GetNFLayer(1);
	for (auto & i : layer2)
	{
		std::cout << i << " ";
		i = 3;
	}
	std::cout << std::endl;

	/*m.SetTBLayer(z - 2, layer);
	m.SetLRLayer(x - 2, layer1);
	m.SetNFLayer(y - 2, layer2);*/

}

int main()
{

#pragma region initial

	using std::cout;
	using std::endl;
	omp_set_num_threads(1);

	// >>>  2D case 

	// !!! Можно создать общую универсальную функцию VonNeumann в которую просто передать номера компонеет котрые надо заполнять наверно с внутренним switchem по знакам между присваиванием!!

	int X{ 70 }; //100
	int Y{ 70 };
	Fluid f(Y, X);
	Medium m(Y, X);

	m.AddCircleInMedium(35, 35, 10);

	f.AddImmersedBodies(m);
	


	//f.Poiseuille_IC(0.01);

	SRTsolver solver(1.0, m, f);
	solver.solve(301);

	// >>>

#pragma endregion

	// >>> 3D case

	//int x{ 10 }; // 10
	//int y{ 35 }; // 35
	//int z{ 10 }; // 10

	////MatrixTest();

	//Fluid3D f(z, y, x);
	//Medium3D m(z,y,x);

	//SRT3DSolver srt(1.0, m, f);
	//srt.solve(51);

	// >>> 

	/*
		- Продумать структуру для BC!
		- Реализовать ГУ Фон-Неймана для всех границ (TOP, BOTTOM, RIGHT)
		- Не распечатывает файл, т.к. нет дирректории DATA
	*/

	return 0;
}
