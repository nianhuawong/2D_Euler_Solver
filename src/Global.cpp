#include "2D_Euler_Solver.h"
#include "QlQr_Solver.h"
#include "Global.h"
#include "Geometry.h"

int num_of_prim_vars;
int current_step, max_num_of_steps;
double cfl_num, time_step;
int method_of_half_q;
int method_of_limiter;
int method_of_flux;
double muscl_k;
double entropy_fix_coeff;
char solve_direction;
int residual_output_steps;
int flow_save_steps;
double converge_criterion;

void Init_Global_Param()
{
	num_of_prim_vars = 4;		//原始变量个数，控制方程个数

	max_num_of_steps = 10000;

	cfl_num   = 0.1;
	time_step = 0.0;			//时间步长要根据最大特征值确定，这里只是初始化

	method_of_half_q  = 1;		//1-MUSCL,		2-WENO,		3-WCNS
	muscl_k			  = 0.0;	//0.0-二阶迎风偏置，		1/3-二阶迎风偏置
	method_of_limiter = 1;		//1-vanleer,	2-minmod,	3-superbee	
	method_of_flux    = 1;		//1-Roe,		2-WENO,		3-WCNS
	entropy_fix_coeff = 0.01;	//Roe格式熵修正系数epsilon

	num_grid_point_x = 101;
	num_grid_point_y = 51;

	solve_direction = 'x';

	residual_output_steps = 20;		//残差输出间隔步数
	flow_save_steps		  = 100;	//流场输出间隔步数
	converge_criterion	  = 1e-8;	//残差收敛标准
}

void Init_Flow()
{
	//流场初始化
	Allocate_3D_Vector(qField,    num_of_prim_vars, total_points_x, total_points_y);
	Allocate_3D_Vector(qField_N1, num_of_prim_vars, total_points_x, total_points_y);
	Allocate_3D_Vector(qField_N2, num_of_prim_vars, total_points_x, total_points_y);
	Allocate_3D_Vector(qField_N3, num_of_prim_vars, total_points_x, total_points_y);	

	//流场赋初值
	vector< vector< int > >& marker = mesh->Get_Marker();

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			if (marker[i][j] == 0) continue;

			qField[IR][i][j] = 1.0;
			qField[IU][i][j] = 3.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = 0.71429;
		}
	}
}

void Compute_Boundary()
{
	vector< vector< int > >& marker = mesh->Get_Marker();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	//左边界：超声速入口
	for (int i = 0; i < ist; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			qField[IR][i][j] = 1.0;
			qField[IU][i][j] = 3.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = 0.71429;
		}
	}

	//上边界：outflow
	for (int i = ist; i < ied; i++)
	{
		for (int j = jed; j < total_points_y; j++)
		{
			qField[IR][i][j] = qField[IR][i][j - 1];
			qField[IU][i][j] = qField[IU][i][j - 1];
			qField[IV][i][j] = qField[IV][i][j - 1];
			qField[IP][i][j] = qField[IP][i][j - 1];
		}
	}

	//右边界：outflow
	for (int i = ied; i < total_points_x; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			if (marker[i - 1][j] == 0) continue;

			qField[IR][i][j] = qField[IR][i - 1][j];
			qField[IU][i][j] = qField[IU][i - 1][j];
			qField[IV][i][j] = qField[IV][i - 1][j];
			qField[IP][i][j] = qField[IP][i - 1][j];
		}
	}

	//下边界, no-slip wall
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst - 1; j >= 0; j--)
		{
			qField[IR][i][j] = qField[IR][i][j + 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j + 1];
		}
	}

	//左物面, no-slip wall
	for (int i = ist + Iw; i < ist + Iw + num_ghost_point; i++)
	{
		for (int j = jst + Jw1; j < jst + Jw2; j++)
		{
			qField[IR][i][j] = qField[IR][i - 1][j];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i - 1][j];
		}
	}

	//上物面, no-slip wall
	for (int i = ist + Iw; i < ied; i++)
	{
		//for (int j = Jw2; j < num_ghost_point + Jw2; j++)
		for (int j = jst + Jw2 - 1; j >= Jw2; j--)
		{
			qField[IR][i][j] = qField[IR][i][j + 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j + 1];
		}
	}

	//下物面, no-slip wall
	for (int i = ist + Iw; i < ied; i++)
	{
		for (int j = jst + Jw1; j < jst + Jw1 + num_ghost_point; j++)
		{
			qField[IR][i][j] = qField[IR][i][j - 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j - 1];
		}
	}

	qField_N1 = qField;
}

void Load_Q()
{
	qField = qField_N1;
}

void Set_Solve_Direction(char direction)
{
	solve_direction = direction;
}

double Energy_2_Pressure(double E, double rho, double u, double v)
{
	double gama = 1.4;
	return (gama - 1) * (E - 0.5 * rho * (u * u + v * v));
}

//计算点的标号范围
void Get_IJK_Region(int& ist, int& ied, int& jst, int& jed)
{
	ist = num_ghost_point;
	ied = num_ghost_point + num_grid_point_x;

	jst = num_ghost_point;
	jed = num_ghost_point + num_grid_point_y;
}