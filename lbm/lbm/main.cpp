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
	Matrix<double> line = m.get_row(2);
	cout << line;
	Matrix<double> row = m.get_coll(2);
	cout << row;
	m.set_coll(3, row);
	cout << m;*/


	// ---------------------------

	
	Fluid f(5, 10);
	Medium m(5, 10);

	f.Poiseuille_IC(0.01);

	SRTsolver solver(1.0, m, f);
	solver.solve();

	return 0;
}
