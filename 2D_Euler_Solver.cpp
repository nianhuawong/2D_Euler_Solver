#include <vector>
#include <iostream>
#include "2D_Euler_Solver.h"
#include "Compute_Half_Node_Q.h"
#include "Compute_Half_Node_Flux.h"
#include "Global_Variables.h"
using namespace std;
using namespace GLOBAL;

int main(int argc, char ** argv )
{
	Simulation * two_dim_Euler_Solver = new Simulation();

	two_dim_Euler_Solver->Run();

	delete two_dim_Euler_Solver;

	return 0;
}

void Simulation::Run()
{
	//建立全局参数、变量、流场
	Init_Global_Param();

	Generate_Mesh();

	Flow_Initialization();

	Compute_Boundary();

	for (current_step = 0; current_step < max_num_of_steps; ++current_step)
	{
		//先计算x方向
		//Load_Q();

		Half_Node_Q();

		Half_Node_Flux();

		//Spatial_Derivative();

		//Time_Integral();
	}

	//再计算y方向

}

