#include "2D_Euler_Solver.h"
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

	Init_Flow();

	for (current_step = 1; current_step <= max_num_of_steps; ++current_step)
	{
		if (current_step == 28)
		{
			int kkk = 1;
		}

		Solve_Time_Step();

		Compute_Boundary();

		Set_Solve_Direction('x');
		Time_Integration();

		Set_Solve_Direction('y');
		Time_Integration();

		Post_Solve();

		if (Need_Stop_Iteration())
		{
			break;
		}

		Set_Field();
	}

	Output_Flowfield();	
}
