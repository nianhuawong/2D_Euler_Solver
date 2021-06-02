#include <cmath>
#include "Time_Step.h"
#include "Geometry.h"
#include "Global.h"
#include "QlQr_Solver.h"

void Solve_Time_Step()
{
	Time_Step* time_step_solver = new Time_Step();
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
#ifdef _OPENMP
#pragma omp parallel for shared( a_max )
#endif
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			if (marker[i][j] == 0) continue;

			double rho = qField[i][j][IR];
			double u   = qField[i][j][IU];
			double v   = qField[i][j][IV];
			double p   = qField[i][j][IP];

			double a = sqrt(fabs(gama * p / rho));
			//a_max = max(a, a_max);
#ifdef _OPENMP
#pragma omp critical
			{
				a_max = max(a, a_max);
			}
#else
			a_max = max(a, a_max);
#endif // _OPENMP	
		}
	}

	double ts1 = cfl_num * dx / a_max;
	double ts2 = cfl_num * dy / a_max;
	time_step = min(ts1, ts2);

	physical_time += time_step;
}


