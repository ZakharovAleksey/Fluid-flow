#include<iostream>
#include<ctime>
#include<cstdlib>
#include<omp.h>

#include"math\2d\my_matrix_2d.h"
#include"phys_values\2d\macroscopic_param_2d.h"
#include"phys_values\2d\distribution_func_2d.h"

#include"modeling_area\medium.h"
#include"modeling_area\fluid.h"
#include"solver\ib_srt.h"
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

	using std::cout;
	using std::endl;
	omp_set_num_threads(1);

#pragma region 2D

#pragma region SRT 

	//int X{ 100 }; //100
	//int Y{ 70 };
	//Fluid f(Y, X);
	//Medium m(Y, X);

	//// Add additional static boundaries
	//m.AddCircleTopFalf(20, 1, 15);
	//m.AddCircleBottomFalf(50, 69, 15);
	////m.AddCircleInMedium(100, 35, 10);
	//f.AddImmersedBodies(m);
	//
	//// Start solution
	//SRTsolver solver(1.0, m, f);
	//solver.Solve(129);

#pragma endregion

#pragma region IB-LBM

	int X{ 102 };
	int Y{ 30 };
	Fluid f(Y, X);
	Medium m(Y, X);

	// Create immersed object
	// тут мы пока не обнавл€ем потожение половины узлов UpdatePosition - убрать!! + сощдать погруженный объ€ект - тромб
	std::unique_ptr<ImmersedBody> body(new ImmersedTromb(102, 30, 32, Point(30, 1), 6));

	// Start solution
	IBSolver s(1.0, f, m, std::move(body));
	s.Solve(1501);

#pragma endregion

#pragma endregion

#pragma region 3D

#pragma region SRT

	//int x{ 10 }; // 10
	//int y{ 35 }; // 35
	//int z{ 10 }; // 10

	////MatrixTest();

	//Fluid3D f(z, y, x);
	//Medium3D m(z,y,x);

	//SRT3DSolver srt(1.0, m, f);
	//srt.Solve(101);

#pragma endregion


#pragma endregion

	return 0;
}
