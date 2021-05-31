﻿#include "2D_Euler_Solver.h"
#include "Global.h"

int main(int argc, char ** argv )
{
	Simulation * two_dim_Euler_Solver = new Simulation();

	two_dim_Euler_Solver->Run();

	delete two_dim_Euler_Solver;

	return 0;
}

void Simulation::Run()
{
	Init_Global_Param();

	Generate_Mesh();

	//Init_Flow();
	Init_Flow_Double_Mach();
	//Compute_Boundary();

	for (current_step = 0; current_step < max_num_of_steps; ++current_step)
	{
		//先计算x方向
		if (current_step == 115)
		{
			int kkk = 1;
		}

		Set_Solve_Direction('x');

		Load_Q();

		//Compute_Boundary();
		Compute_Boundary_Double_Mach();

		Time_Integration();

		//再计算y方向
		Set_Solve_Direction('y');

		Time_Integration();

		Compute_Residual();

		Output_Flowfield();

		//Load_Q();

		//Compute_Boundary();

		if (stop_by_residual)
		{
			break;
		}
	}

	Output_Flowfield();	
}

void Test()
{
	vector < vector< double > > A = { {1,2,3,0.5},{3,2,1,1.2},{1,3,2,2.1} };
	//vector < vector< double > > B = { {2,3,4},{3,4,2},{4,2,3},{2,4,7} };
	//vector < vector< double > > B = { {2,3},{3,4},{4,2},{2,4} };

	//vector < vector< double > > C;
	//Allocate_2D_Vector(C,3,2);
	//MatrixMultiply(A,B,C,3,4,2);

	vector< double > B = { 2, 3, 9, 12 };
	vector< double > C(3);
	MatrixMultiply(A, B, C, 3, 4);
}
