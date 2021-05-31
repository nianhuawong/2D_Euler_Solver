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
	Get_IJK_Region(ist, ied, jst, jed);
	Allocate_3D_Vector(rhs, num_half_point_x, num_half_point_y, num_of_prim_vars);
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
	VInt2D& marker = mesh->Get_Marker();

	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;

			VDouble rhsVector = rhs[i][j];
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = -(fluxVector[i][j][iVar] - fluxVector[i - 1][j][iVar]) / dx;
			}
			rhs[i][j] = rhsVector;
		}
	}
}

void Spatial_Derivative::Spatial_Derivative_Y()
{
	VInt2D& marker = mesh->Get_Marker();

	for (int i = ist; i < ied - 1; i++)
	{
		for (int j = jst; j < jed - 1; j++)
		{
			if (marker[i][j] == 0) continue;

			VDouble rhsVector = rhs[i][j];
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = -(fluxVector[i][j][iVar] - fluxVector[i][j - 1][iVar]) / dy;
			}
			rhs[i][j] = rhsVector;
		}
	}
}
