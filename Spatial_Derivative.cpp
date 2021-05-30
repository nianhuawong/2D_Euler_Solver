#include "2D_Euler_Solver.h"
#include "Spatial_Derivative.h"
#include "Geometry.h"
#include "Global.h"
#include "Flux_Solver.h"

vector< vector< vector< double > > > rhs;
void Spatial_Derivative()
{
	rhs.resize(num_of_prim_vars);

	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		Allocate_2D_Vector(rhs[iVar], num_half_point_x, num_half_point_y);
		Allocate_2D_Vector(rhs[iVar], num_half_point_x, num_half_point_y);
	}

	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		for (int j = 0; j < num_half_point_y; j++)
		{
			for (int i = 0; i < num_half_point_x; i++)
			{
				rhs[iVar][i][j] = -(fluxVector[i][j][iVar] - fluxVector[i - 1][j][iVar]) / dx;
			}
		}
	}
}