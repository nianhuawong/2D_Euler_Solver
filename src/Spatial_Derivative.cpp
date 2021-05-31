#include "2D_Euler_Solver.h"
#include "Spatial_Derivative.h"
#include "Geometry.h"
#include "Global.h"
#include "Flux_Solver.h"
#include <iostream>

VDouble3D rhs;
void Solve_Spatial_Derivative()
{
	auto* spatial_derivative = new Spatial_Derivative();

	spatial_derivative->Compute_Spatial_Derivative();

	delete spatial_derivative;
}

Spatial_Derivative::Spatial_Derivative()
{
	rhs.resize(num_of_prim_vars);

	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		Allocate_2D_Vector(rhs[iVar], num_half_point_x, num_half_point_y);		
	}
}

void Spatial_Derivative::Compute_Spatial_Derivative()
{
	if (solve_direction == 'x')
	{
		this->Spatial_Derivative_X();
	}
	else if (solve_direction == 'y')
	{
		this->Spatial_Derivative_Y();
	}
	else
	{
		cout << "³ö´í£¬Çë¼ì²é£¡" << endl;
	}
}

void Spatial_Derivative::Spatial_Derivative_X()
{
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

void Spatial_Derivative::Spatial_Derivative_Y()
{
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		for (int i = 0; i < num_half_point_x; i++)
		{
			for (int j = 0; j < num_half_point_y; j++)
			{
				rhs[iVar][i][j] = -(fluxVector[i][j][iVar] - fluxVector[i - 1][j][iVar]) / dy;
			}
		}
	}
}
