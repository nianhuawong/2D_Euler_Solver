#include <cmath>
#include "2D_Euler_Solver.h"
#include "Geometry.h"
#include "QlQr_Solver.h"

void Compute_Boundary()
{
	if (global_case_id == 1)
	{
		Compute_Boundary_Double_Mach();
	}
	else if (global_case_id == 2)
	{
		Compute_Boundary_Blunt_Body();
	}
}

void Compute_Boundary_Blunt_Body()
{
	VInt2D& marker = mesh->Get_Marker_Q();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	//上边界：outflow
	for (int i = ist; i <= ied; i++)
	{
		for (int j = jed; j <= jed + 2; j++)
		{
			qField[i][j][IR] = qField[i][j - 1][IR];
			qField[i][j][IU] = qField[i][j - 1][IU];
			qField[i][j][IV] = qField[i][j - 1][IV];
			qField[i][j][IP] = qField[i][j - 1][IP];
		}
	}

	//右边界：outflow
	for (int j = jst; j <= jed; j++)
	{
		for (int i = ied; i <= ied + 2; i++)
		{
			if (marker[i][j] == 0) continue;

			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = qField[i - 1][j][IU];
			qField[i][j][IV] = qField[i - 1][j][IV];
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}
	//左边界：超声速入口
	for (int i = 0; i <= ist; i++)
	{
		for (int j = jst; j <= jed; j++)
		{
			qField[i][j][IR] = 1.0;
			qField[i][j][IU] = 3.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = 0.71429;
		}
	}

	//下边界, slip wall
	for (int i = ist; i <= ied; i++)
	{
		qField[i][jst][IR] =  qField[i][jst + 1][IR];		//weno格式不加这一行在下边界不会反射激波？
		//qField[i][jst][IU] =  qField[i][jst + 1][IU];
		qField[i][jst][IV] =  0.0;
		qField[i][jst][IP] =  qField[i][jst + 1][IP];		//weno格式不加这一行在下边界不会反射激波？

		qField[i][jst - 1][IR] =  qField[i][jst + 1][IR];
		qField[i][jst - 1][IU] =  qField[i][jst + 1][IU];
		qField[i][jst - 1][IV] = -qField[i][jst + 1][IV];
		qField[i][jst - 1][IP] =  qField[i][jst + 1][IP];

		qField[i][jst - 2][IR] =  qField[i][jst + 2][IR];
		qField[i][jst - 2][IU] =  qField[i][jst + 2][IU];
		qField[i][jst - 2][IV] = -qField[i][jst + 2][IV];
		qField[i][jst - 2][IP] =  qField[i][jst + 2][IP];
	}

	//上物面, slip wall
	int JJ = jst + Jw2;
	for (int i = ist + Iw; i <= ied; i++)
	{
		//qField[i][JJ][IR] = qField[i][JJ + 1][IR];
		//qField[i][JJ][IU] = qField[i][JJ + 1][IU];
		qField[i][JJ][IV] = 0.0;
		//qField[i][JJ][IP] = qField[i][JJ + 1][IP];

		qField[i][JJ - 1][IR] =  qField[i][JJ + 1][IR];
		qField[i][JJ - 1][IU] =  qField[i][JJ + 1][IU];
		qField[i][JJ - 1][IV] = -qField[i][JJ + 1][IV];
		qField[i][JJ - 1][IP] =  qField[i][JJ + 1][IP];	
	
		qField[i][JJ - 2][IR] =  qField[i][JJ + 2][IR];
		qField[i][JJ - 2][IU] =  qField[i][JJ + 2][IU];
		qField[i][JJ - 2][IV] = -qField[i][JJ + 2][IV];
		qField[i][JJ - 2][IP] =  qField[i][JJ + 2][IP];
	}

	//下物面, slip wall
	JJ = jst + Jw1;
	for (int i = ist + Iw; i <= ied; i++)
	{
		//qField[i][JJ][IR] =  qField[i][JJ - 1][IR];
		//qField[i][JJ][IU] =  qField[i][JJ - 1][IU];
		qField[i][JJ][IV] =  0.0;
		//qField[i][JJ][IP] =  qField[i][JJ - 1][IP];	

		qField[i][JJ + 1][IR] =  qField[i][JJ - 1][IR];
		qField[i][JJ + 1][IU] =  qField[i][JJ - 1][IU];
		qField[i][JJ + 1][IV] = -qField[i][JJ - 1][IV];
		qField[i][JJ + 1][IP] =  qField[i][JJ - 1][IP];	

		qField[i][JJ + 2][IR] =  qField[i][JJ - 2][IR];
		qField[i][JJ + 2][IU] =  qField[i][JJ - 2][IU];
		qField[i][JJ + 2][IV] = -qField[i][JJ - 2][IV];
		qField[i][JJ + 2][IP] =  qField[i][JJ - 2][IP];
	}

	//左物面, slip wall
	int II = ist + Iw;
	for (int j = jst + Jw1; j <= jst + Jw2; j++)
	{
		//qField[II][j][IR] =  qField[II - 1][j][IR];
		qField[II][j][IU] = 0.0;
		//qField[II][j][IV] =  0.0; //qField[II - 1][j][IV];
		//qField[II][j][IP] =  qField[II - 1][j][IP];

		qField[II + 1][j][IR] =  qField[II - 1][j][IR];
		qField[II + 1][j][IU] = -qField[II - 1][j][IU];
		qField[II + 1][j][IV] =  qField[II - 1][j][IV];
		qField[II + 1][j][IP] =  qField[II - 1][j][IP];

		qField[II + 2][j][IR] =  qField[II - 2][j][IR];
		qField[II + 2][j][IU] = -qField[II - 2][j][IU];
		qField[II + 2][j][IV] =  qField[II - 2][j][IV];
		qField[II + 2][j][IP] =  qField[II - 2][j][IP];
	}

	//角点处理，左上，设置参考角点标号
	II = ist + Iw;
	JJ = jst + Jw2;
	//1点，对角方向
	qField[II + 2][JJ - 2][IR] =  qField[II][JJ][IR];
	qField[II + 2][JJ - 2][IU] = -qField[II][JJ][IU];
	qField[II + 2][JJ - 2][IV] = -qField[II][JJ][IV];
	qField[II + 2][JJ - 2][IP] =  qField[II][JJ][IP];

	//2点, x方向取值
	qField[II + 1][JJ - 2][IR] =  qField[II - 1][JJ - 2][IR];
	qField[II + 1][JJ - 2][IU] = -qField[II - 1][JJ - 2][IU];
	qField[II + 1][JJ - 2][IV] =  qField[II - 1][JJ - 2][IV];
	qField[II + 1][JJ - 2][IP] =  qField[II - 1][JJ - 2][IP];

	//3点, y方向取值
	qField[II + 2][JJ - 1][IR] =  qField[II + 2][JJ + 1][IR];
	qField[II + 2][JJ - 1][IU] =  qField[II + 2][JJ + 1][IU];
	qField[II + 2][JJ - 1][IV] = -qField[II + 2][JJ + 1][IV];
	qField[II + 2][JJ - 1][IP] =  qField[II + 2][JJ + 1][IP];

	//4点， x方向取值和y方向取值的算术平均
	qField[II + 1][JJ - 1][IR] = 0.5 * (qField[II - 1][JJ - 1][IR] + qField[II + 1][JJ + 1][IR]);
	qField[II + 1][JJ - 1][IU] = 0.5 * (qField[II - 1][JJ - 1][IU] + qField[II + 1][JJ + 1][IU]);
	qField[II + 1][JJ - 1][IV] = 0.5 * (qField[II - 1][JJ - 1][IV] + qField[II + 1][JJ + 1][IV]);
	qField[II + 1][JJ - 1][IP] = 0.5 * (qField[II - 1][JJ - 1][IP] + qField[II + 1][JJ + 1][IP]);
	
	//=========================
	//角点处理，左下，设置参考角点标号
	II = ist + Iw;
	JJ = jst + Jw1;

	//1点，对角方向
	qField[II + 2][JJ + 2][IR] =  qField[II][JJ][IR];
	qField[II + 2][JJ + 2][IU] = -qField[II][JJ][IU];
	qField[II + 2][JJ + 2][IV] = -qField[II][JJ][IV];
	qField[II + 2][JJ + 2][IP] =  qField[II][JJ][IP];

	//2点, x方向取值
	qField[II + 1][JJ + 2][IR] =  qField[II - 1][JJ + 2][IR];
	qField[II + 1][JJ + 2][IU] = -qField[II - 1][JJ + 2][IU];
	qField[II + 1][JJ + 2][IV] =  qField[II - 1][JJ + 2][IV];
	qField[II + 1][JJ + 2][IP] =  qField[II - 1][JJ + 2][IP];

	//3点, y方向取值
	qField[II + 2][JJ + 1][IR] =  qField[II + 2][JJ - 1][IR];
	qField[II + 2][JJ + 1][IU] =  qField[II + 2][JJ - 1][IU];
	qField[II + 2][JJ + 1][IV] = -qField[II + 2][JJ - 1][IV];
	qField[II + 2][JJ + 1][IP] =  qField[II + 2][JJ - 1][IP];

	//4点， x方向取值和y方向取值的算术平均
	qField[II + 1][JJ + 1][IR] = 0.5 * (qField[II - 1][JJ + 1][IR] + qField[II + 1][JJ - 1][IR]);
	qField[II + 1][JJ + 1][IU] = 0.5 * (qField[II - 1][JJ + 1][IU] + qField[II + 1][JJ - 1][IU]);
	qField[II + 1][JJ + 1][IV] = 0.5 * (qField[II - 1][JJ + 1][IV] + qField[II + 1][JJ - 1][IV]);
	qField[II + 1][JJ + 1][IP] = 0.5 * (qField[II - 1][JJ + 1][IP] + qField[II + 1][JJ - 1][IP]);
}

void Compute_Boundary_Double_Mach()
{
	VInt2D& marker = mesh->Get_Marker_Q();
	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	double ssw = 10.0 / sin(PI / 3);
	double xsw = 1.0 / 6 + 1.0 / tan(PI / 3) + ssw * physical_time;

	//左边界：inflow
	//for (int i = 0; i < ist; i++)
	for (int i = 0; i <= ist; i++)
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
	//for (int j = jed + 1; j <= jed + 2; j++)
	for (int j = jed; j <= jed + 2; j++)
	{
		for (int i = ist; i <= ied; i++)
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
	for (int j = jst; j <= jed; j++)
	{
		//for (int i = ied + 1; i <= ied + 2; i++)
		for (int i = ied; i <= ied + 2; i++)
		{
			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = qField[i - 1][j][IU];
			qField[i][j][IV] = qField[i - 1][j][IV];
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}

	//下边界
	for (int i = ist; i <= ied; i++)
	{
		//for (int j = jst - 1; j >= 0; j--)
		for (int j = jst; j >= 0; j--)
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
