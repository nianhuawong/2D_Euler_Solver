// 2D_Euler_Solver.cpp: 定义应用程序的入口点。
//

#include "2D_Euler_Solver.h"

using namespace std;

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


	for (iter = 0; iter < numberOfTimeSteps; ++iter)
	{
		//先计算x方向
		//Load_Q();

		//Half_Node_Q();

		//Half_Node_Flux();

		//Spatial_Derivative();

		//Time_Integral();
	}

	//再计算y方向

}