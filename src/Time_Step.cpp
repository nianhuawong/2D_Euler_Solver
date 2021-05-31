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
	Get_IJK_Region(ist, ied, jst, jed);
}

void Time_Step::Compute_Time_Step()
{
	double a_max = -1e40;
	double gama  = 1.4;

	VInt2D& marker = mesh->Get_Marker();
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			if (marker[i][j] == 0) continue;

			double rho = qField[IR][i][j];
			double u = qField[IU][i][j];
			double v = qField[IV][i][j];
			double p = qField[IP][i][j];

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


