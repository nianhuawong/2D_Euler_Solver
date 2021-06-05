#include "2D_Euler_Solver.h"
#include "Spatial_Derivative.h"
#include "Geometry.h"
#include "Global.h"
#include "Flux_Solver.h"
#include <iostream>

VDouble3D rhs;
void Solve_Spatial_Derivative()
{
	Spatial_Derivative* spatial_derivative = new Spatial_Derivative();

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
	if (method_of_flux == 4)
	{
		if (solve_direction == 'x')
		{
			this->Spatial_Derivative_WCNS_X();
		}
		else if (solve_direction == 'y')
		{
			this->Spatial_Derivative_WCNS_Y();
		}
	}
	else
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
			cout << "出错，请检查！" << endl;
		}
	}

}

void Spatial_Derivative::Spatial_Derivative_X()
{
	VInt2D& marker = mesh->Get_Marker();

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 1; i <= ied + 1; i++)//i=0和i=ied+2的空间导数项没有值
	{
		for (int j = jst; j <= jed; j++)
		{
			if (marker[i][j] == 0) continue;
			if (i == 61 && j == 62)
			{
				int kkk = 1;
			}
			VDouble rhsVector(num_of_prim_vars);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = -(fluxVector[i][j][iVar] - fluxVector[i - 1][j][iVar]) / dx;
			}
			rhs[i][j] = rhsVector;
			if (IsNaN(rhsVector))
			{
				int kkk = 1;
			}
		}
	}
}

void Spatial_Derivative::Spatial_Derivative_Y()
{
	VInt2D& marker = mesh->Get_Marker();

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 1; j <= jed + 1; j++)//j=0和j=jed+2的空间导数项没有值
		{
			if (marker[i][j] == 0) continue;
			if (i == 16 && j == 2)//((i==60||i == 61) && j == 62)
			{
				int kkk = 1;
			}
			VDouble rhsVector(num_of_prim_vars);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = -(fluxVector[i][j][iVar] - fluxVector[i][j - 1][iVar]) / dy;
			}
			rhs[i][j] = rhsVector;
			if (IsNaN(rhsVector))
			{
				int kkk = 1;
			}
		}
	}
}

void Spatial_Derivative::Spatial_Derivative_WCNS_X()
{
	double ds = dx;
	double a  = 75.0 / 64.0, b = -25.0 / 384.0, c = 3.0 / 640;

	VInt2D& marker = mesh->Get_Marker();
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = ist + 1; i <= ied - 1; i++)	//i=0, 1, 2，ied, ied+1, ied+2没有算
		{
			if (marker[i][j] == 0) continue;

			VDouble rhsVector(num_of_prim_vars);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = a / ds * (fluxVector[i    ][j][iVar] - fluxVector[i - 1][j][iVar])
								+ b / ds * (fluxVector[i + 1][j][iVar] - fluxVector[i - 2][j][iVar])
								+ c / ds * (fluxVector[i + 2][j][iVar] - fluxVector[i - 3][j][iVar]);

				rhsVector[iVar] = - rhsVector[iVar];
			}
			rhs[i][j] = rhsVector;

			if (IsNaN(rhsVector))
			{
				int kkk = 1;
			}
		}
	}
	//边界处导数i=1, 2，ied, ied+1.
	for (int j = jst; j <= jed; j++)
	{
		for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
		{
			rhs[1][j][iVar] =-(-11.0/12 * fluxVector[0][j][iVar] + 17.0/24 * fluxVector[1][j][iVar]
							   + 3.0/8  * fluxVector[2][j][iVar] -  5.0/24 * fluxVector[3][j][iVar]
							   + 1.0/24 * fluxVector[4][j][iVar]) / ds;

			rhs[2][j][iVar] =-(- 1.0/24 * fluxVector[0][j][iVar] -  9.0/8  * fluxVector[1][j][iVar]
							   + 9.0/8	* fluxVector[2][j][iVar] -  1.0/24 * fluxVector[3][j][iVar]) / ds;

			rhs[ied    ][j][iVar] =-(- 1.0/24 * fluxVector[ied + 1][j][iVar] -  9.0/8  * fluxVector[ied    ][j][iVar]
									 + 9.0/8  * fluxVector[ied - 1][j][iVar] -  1.0/24 * fluxVector[ied - 2][j][iVar]) / ds;

			rhs[ied + 1][j][iVar] =-(-11.0/12 * fluxVector[ied + 1][j][iVar] + 17.0/24 * fluxVector[ied    ][j][iVar]
									 + 3.0/8  * fluxVector[ied - 1][j][iVar] -  5.0/24 * fluxVector[ied - 2][j][iVar]
									 + 1.0/24 * fluxVector[ied - 3][j][iVar]) / ds;
		}
	}
}

void Spatial_Derivative::Spatial_Derivative_WCNS_Y()
{
	double ds = dy;
	double a  = 75.0 / 64.0, b = -25.0 / 384.0, c = 3.0 / 640;
	
	VInt2D& marker = mesh->Get_Marker();
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = jst + 1; j <= jed - 1; j++)//j=0, 1, 2，jed, jed+1, jed+2没有算
		{
			if (marker[i][j] == 0) continue;

			VDouble rhsVector(num_of_prim_vars);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = a / ds * (fluxVector[i][j    ][iVar] - fluxVector[i][j - 1][iVar])
								+ b / ds * (fluxVector[i][j + 1][iVar] - fluxVector[i][j - 2][iVar])
								+ c / ds * (fluxVector[i][j + 2][iVar] - fluxVector[i][j - 3][iVar]);

				rhsVector[iVar] = - rhsVector[iVar];
			}
			rhs[i][j] = rhsVector;
		}
	}
	//边界处导数j=1, 2，jed, jed+1
	for (int i = ist; i <= ied; i++)
	{
		for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
		{
			rhs[i][1][iVar] =-(-11.0/12 * fluxVector[i][0][iVar] + 17.0/24 * fluxVector[i][1][iVar]
							   + 3.0/8  * fluxVector[i][2][iVar] -  5.0/24 * fluxVector[i][3][iVar]
							   + 1.0/24 * fluxVector[i][4][iVar]) / ds;

			rhs[i][2][iVar] =-(- 1.0/24 * fluxVector[i][0][iVar] -  9.0/8  * fluxVector[i][1][iVar]
							   + 9.0/8	* fluxVector[i][2][iVar] -  1.0/24 * fluxVector[i][3][iVar]) / ds;

			rhs[i][jed    ][iVar] =-(- 1.0/24 * fluxVector[i][jed + 1][iVar] -  9.0/8  * fluxVector[i][jed    ][iVar]
									 + 9.0/8  * fluxVector[i][jed - 1][iVar] -  1.0/24 * fluxVector[i][jed - 2][iVar]) / ds;

			rhs[i][jed + 1][iVar] =-(-11.0/12 * fluxVector[i][jed + 1][iVar] + 17.0/24 * fluxVector[i][jed    ][iVar]
									 + 3.0/8  * fluxVector[i][jed - 1][iVar] -  5.0/24 * fluxVector[i][jed - 2][iVar]
									 + 1.0/24 * fluxVector[i][jed - 3][iVar]) / ds;
		}
	}
}
