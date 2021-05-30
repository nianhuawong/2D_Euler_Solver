#include "Time_Step.h"
#include "Geometry.h"
#include "Global.h"
#include "QlQr_Solver.h"

void Solve_Time_Step()
{
	auto* time_step_solver = new Time_Step();
	time_step_solver->Compute_Time_Step();
	delete time_step_solver;
}

Time_Step::Time_Step()
{

}

void Time_Step::Compute_Time_Step()
{
	double a_max = -1e40;
	double gama  = 1.4;
	for (int i = 0; i < num_half_point_x; i++)
	{
		for (int j = 0; j < num_half_point_y; j++)
		{
			double rho = qField1[IR][i][j];
			double u = qField1[IU][i][j];
			double v = qField1[IV][i][j];
			double E = qField1[IP][i][j];
			double p = Energy_2_Pressure(E, rho, u, v);

			double a = gama * p / rho;

			a_max = a > a_max ? a : a_max;
		}
	}

	if (solve_direction == 'x')
	{
		time_step = cfl_num * dx / a_max;
	}
	else if (solve_direction == 'y')
	{
		time_step = cfl_num * dy / a_max;
	}
}


