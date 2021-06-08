#include <cmath>
#include "2D_Euler_Solver.h"
#include "Geometry.h"
#include "QlQr_Solver.h"

void Init_Flow()
{
	if (global_case_id == 1)
	{
		Init_Flow_Double_Mach();
	}
	else if (global_case_id == 2)
	{
		Init_Flow_Blunt_Body();
	}
}

void Init_Flow_Blunt_Body()
{
	//流场初始化
	Allocate_3D_Vector(qField, total_points_x, total_points_y, num_of_prim_vars);
	Allocate_3D_Vector(qField_N0, total_points_x, total_points_y, num_of_prim_vars);
	Allocate_3D_Vector(qField_N1, total_points_x, total_points_y, num_of_prim_vars);

	//流场赋初值
	VInt2D& marker = mesh->Get_Marker_Q();

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	for (int i = ist; i <= ied; i++)
	{
		for (int j = jst; j <= jed; j++)
		{
			if (marker[i][j] == 0) continue;

			qField[i][j][IR] = 1.0;
			qField[i][j][IU] = 3.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = 0.71429;
		}
	}
}

void Init_Flow_Double_Mach()
{
	//流场初始化
	Allocate_3D_Vector(qField, total_points_x, total_points_y, num_of_prim_vars);
	Allocate_3D_Vector(qField_N1, total_points_x, total_points_y, num_of_prim_vars);

	//流场赋初值
	VInt2D& marker = mesh->Get_Marker_Q();
	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();

	double gama = 1.4;
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	for (int i = ist; i <= ied; i++)
	{
		for (int j = jst; j <= jed; j++)
		{
			double x_node, y_node;

			grid_points[i][j].Get_Point_Coord(x_node, y_node);
			if (x_node <= (1.0 / 6.0 + y_node / tan(PI / 3)))
			{
				qField[i][j][IR] = 8.0;
				qField[i][j][IU] = 8.25 * cos(PI / 6);
				qField[i][j][IV] = -8.25 * sin(PI / 6);
				qField[i][j][IP] = 116.5;
			}
			else if (x_node > (1.0 / 6.0 + y_node / tan(PI / 3)))
			{
				qField[i][j][IR] = 1.4;
				qField[i][j][IU] = 0.0;
				qField[i][j][IV] = 0.0;
				qField[i][j][IP] = 1.0;
			}
		}
	}
}
