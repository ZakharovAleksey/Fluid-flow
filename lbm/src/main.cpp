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
#include"solver\mrt.h"
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

	int X{ 30 }; //100
	int Y{ 32 };
	Medium m(Y, X);
	//m.AddCircle(15, 15, 6);
	Fluid f(Y, X, m);
	BCs bc(f.f_);

	bc.SetBounceBackBC(Boundary::TOP);
	bc.SetBounceBackBC(Boundary::BOTTOM);
	bc.SetPeriodicBC(Boundary::LEFT, Boundary::RIGHT);

	// Set vx parabolic profile
	std::vector<double> vx;
	for (int i = 1; i < Y - 1; ++i)
	{
		double val = 0.01; // *sin((i - 1) * M_PI / (Y - 3));
		vx.push_back(val);
	}
	vx.at(vx.size() - 1) = 0.0;


	//bc.SetVonNeumannBC(Boundary::LEFT, vx, std::vector<double>(vx.size(), 0.0));
	//bc.SetVonNeumannBC(Boundary::RIGHT, std::vector<double>(vx.size(), 0.0), std::vector<double>(vx.size(), 0.0));

	//bc.SetDirichletBC(Boundary::LEFT, std::vector<double>(Y - 2, 1.001));
	//bc.SetDirichletBC(Boundary::RIGHT, std::vector<double>(Y -2, 1.0));
	
	//// Start solution
	//SRTsolver solver(1.0, m, f, &bc);
	//solver.Solve(15001);

#pragma endregion

#pragma region MRT

	//int X{ 200 }; //100
	//int Y{ 40 };
	//Medium m(Y, X);
	//Fluid f(Y, X, m);

#pragma region obstacles

	///* Y-shape channel flow obstacles
	//m.AddSquare(1, Y-2, X/2, Y/4);
	//m.AddTopAngle(X / 2, 1, Y / 4);
	//m.AddSquare(1, Y/4, X/2, Y/4);
	//m.AddBottomAngle(X / 2, Y - 2, Y / 4);
	//m.AddSquare(X / 2 + Y / 2, 3 * Y / 4, X / 2, Y / 2);
	//m.AddCircle(X / 2 + Y / 2, Y / 2, Y / 4);*/
	// >>

	/* >> Around tromb flow obstacle
	m.AddCircleTopFalf(35, 39, 15);
	 >>*/
	
	/* >> Flow around cylinder obstacle
	m.AddCircle(X/10, Y/2, Y/10+ 1);
	 >>*/

	/* >> Flow in asterios obstacle
	m.AddSquare(1, 0.25 * Y, X / 2 - Y / 4, 0.25 * Y);
	m.AddTopAngle(X / 2 - Y / 4, -1, Y / 4);
	m.AddRightAngle(X / 2, -1, Y / 4);
	m.AddSquare(X / 2 + Y /4, 0.25 * Y, X / 2 - Y / 4, 0.25 * Y);
	>>*/

	 /*>> Flow between trombs obstacles
	 m.AddCircle(35, 39, 10);
	 m.AddCircle(35, 1, 10);
	 >> */

#pragma endregion


	// Start solution
	//MRTSolver solver(0.514, m, f, &bc); // tau = 0.5008
	//solver.Solve(100001);

#pragma endregion


#pragma region IB-LBM
	
	ImmersedBody* rbc(new ImmersedRBC(X, Y, 32, Point(15, 15), 6));
	//ImmersedBody* tromb(new ImmersedCircle(X, Y, 64, Point(38, 35), 15, M_PI, 2.0 * M_PI));
	std::vector<ImmersedBody*> bodies;
	//bodies.push_back(rbc);

	// Start solution
	IBSolver s(1.0, f, m, &bc, bodies);
	s.Solve(2001);


	//IBMRTSolver sol(0.6, f, m, &bc, bodies);
	//sol.Solve(1501);

#pragma endregion

#pragma endregion

#pragma region 3D

#pragma region SRT

	//int x{ 12 }; // 10
	//int y{ 12 }; // 35
	//int z{ 32 }; // 10

	////MatrixTest();

	//Fluid3D f(z, y, x);
	//Medium3D m(z,y,x);

	//SRT3DSolver srt(1.0, m, f);
	//srt.Solve(51);

#pragma endregion


#pragma endregion

	return 0;
}
