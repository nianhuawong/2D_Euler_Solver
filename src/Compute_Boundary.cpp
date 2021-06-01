#include "2D_Euler_Solver.h"
#include "Geometry.h"
#include "QlQr_Solver.h"

void Compute_Boundary()
{
	//Compute_Boundary_Blunt_Body();
	Compute_Boundary_Double_Mach();
}

void Compute_Boundary_Blunt_Body()
{
	VInt2D& marker = mesh->Get_Marker();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	//左边界：超声速入口
	for (int i = 0; i < ist; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			qField[i][j][IR] = 1.0;
			qField[i][j][IU] = 3.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = 0.71429;
		}
	}

	//上边界：outflow
	for (int i = ist; i < ied; i++)
	{
		for (int j = jed; j < total_points_y; j++)
		{
			qField[i][j][IR] = qField[i][j - 1][IR];
			qField[i][j][IU] = qField[i][j - 1][IU];
			qField[i][j][IV] = qField[i][j - 1][IV];
			qField[i][j][IP] = qField[i][j - 1][IP];
		}
	}

	//右边界：outflow
	for (int j = jst; j < jed; j++)
	{
		for (int i = ied; i < total_points_x; i++)
		{
			if (marker[i - 1][j] == 0) continue;

			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = qField[i - 1][j][IU];
			qField[i][j][IV] = qField[i - 1][j][IV];
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}

	//下边界, no-slip wall
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst - 1; j >= 0; j--)
		{
			qField[i][j][IR] = qField[i][j + 1][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i][j + 1][IP];
		}
	}

	//左物面, no-slip wall
	for (int j = jst + Jw1; j < jst + Jw2; j++)
	{
		for (int i = ist + Iw; i < ist + Iw + num_ghost_point; i++)
		{
			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}

	//上物面, no-slip wall
	for (int i = ist + Iw; i < ied; i++)
	{
		//for (int j = Jw2; j < num_ghost_point + Jw2; j++)
		for (int j = jst + Jw2 - 1; j >= Jw2; j--)
		{
			qField[i][j][IR] = qField[i][j + 1][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i][j + 1][IP];
		}
	}

	//下物面, no-slip wall
	for (int i = ist + Iw; i < ied; i++)
	{
		for (int j = jst + Jw1; j < jst + Jw1 + num_ghost_point; j++)
		{
			qField[i][j][IR] = qField[i][j - 1][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i][j - 1][IP];
		}
	}
}

void Compute_Boundary_Double_Mach()
{
	VInt2D& marker = mesh->Get_Marker();
	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	double ssw = 10.0 / sin(PI / 3);
	double xsw = 1.0 / 6 + 1.0 / tan(PI / 3) + ssw * physical_time;

	//左边界：inflow
	for (int i = 0; i < ist; i++)
	//for (int i = 0; i <= ist; i++)
	{
		for (int j = jst; j <= jed; j++)
		{
			qField[i][j][IR] = 8.0;
			qField[i][j][IU] = 8.25 * cos(PI / 6);
			qField[i][j][IV] = -8.25 * sin(PI / 6);
			qField[i][j][IP] = 116.5;
		}
	}

	//上边界
	for (int j = jed; j < total_points_y; j++)
	//for (int j = jed - 1; j < total_points_y; j++)
	{
		for (int i = ist; i < ied; i++)
		{
			double x_node, y_node;
			grid_points[i][j].Get_Point_Coord(x_node, y_node);
			if (x_node >= 0 && x_node <= xsw)
			{
				qField[i][j][IR] = 8.0;
				qField[i][j][IU] = 8.25 * cos(PI / 6);
				qField[i][j][IV] = -8.25 * sin(PI / 6);
				qField[i][j][IP] = 116.5;
			}
			else if (x_node > xsw && x_node <= 4)
			{
				qField[i][j][IR] = 1.4;
				qField[i][j][IU] = 0.0;
				qField[i][j][IV] = 0.0;
				qField[i][j][IP] = 1.0;
			}
		}
	}

	//右边界：outflow
	for (int j = jst; j < jed; j++)
	{
		for (int i = ied; i < total_points_x; i++)
		//for (int i = ied - 1; i < total_points_x; i++)
		{
			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = qField[i - 1][j][IU];
			qField[i][j][IV] = qField[i - 1][j][IV];
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}

	//下边界
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst - 1; j >= 0; j--)
		//for (int j = jst; j >= 0; j--)
		{
			double x_node, y_node;
			grid_points[i][j].Get_Point_Coord(x_node, y_node);
			if (x_node >= 0 && x_node <= 1.0 / 6)
			{
				qField[i][j][IR] = 8.0;
				qField[i][j][IU] = 8.25 * cos(PI / 6);
				qField[i][j][IV] = -8.25 * sin(PI / 6);
				qField[i][j][IP] = 116.5;
			}
			else if (x_node > 1.0 / 6 && x_node <= 4)
			{
				qField[i][j][IR] = qField[i][j + 1][IR];
				qField[i][j][IU] = qField[i][j + 1][IU];
				qField[i][j][IV] = 0.0;
				qField[i][j][IP] = qField[i][j + 1][IP];
			}
		}
	}
}
