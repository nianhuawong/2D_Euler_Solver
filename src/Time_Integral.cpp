#include "2D_Euler_Solver.h"
#include "Spatial_Derivative.h"
#include "Geometry.h"
#include "Global.h"
#include "Flux_Solver.h"
#include "QlQr_Solver.h"
#include "Time_Integral.h"


void Time_Integration()
{
	auto* time_marching = new Time_Marching_Solver();
	time_marching->Time_Marching();
	delete time_marching;
}

Time_Marching_Solver::Time_Marching_Solver()
{

}

void Time_Marching_Solver::Time_Marching()
{
	Load_Q();

	Solve_QlQr();

	Solve_Flux();

	Solve_Spatial_Derivative();

	Solve_Time_Step();

	Update_Flowfield();
}

void Update_Flowfield()
{
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		for (int j = 0; j < num_half_point_y; j++)
		{
			for (int i = 0; i < num_half_point_x; i++)
			{
				qField_N1[iVar][i][j] = qField[iVar][i][j] + time_step * rhs[iVar][i][j];
			}
		}
	}

//qField_N2[iVar][i][j] = 3.0 / 4.0 * qField[iVar][i][j] + 1.0 / 4.0 * qField_N1[iVar][i][j] + 1.0 / 4.0 * time_step * rhs1[iNode];

//qField_N3[iVar][i][j] = 1.0 / 3.0 * qField[iVar][i][j] + 2.0 / 3.0 * qField_N2[iVar][i][j] + 2.0 / 3.0 * time_step * rhs2[iNode];

}