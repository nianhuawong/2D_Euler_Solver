#include "2D_Euler_Solver.h"
#include "Compute_Half_Node_Q.h"
#include "Global_Variables.h"
#include "Geometry.h"
using namespace GLOBAL;

void GLOBAL::Init_Global_Param()
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
}

void GLOBAL::Flow_Initialization()
{
	//流场初始化
	qField.resize(num_of_prim_vars);
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		Allocate_2D_Vector(qField[iVar], total_points_x, total_points_y);
	}

	//流场赋初值
	vector< vector< int > >& marker = mesh->Get_Marker();
	for (int i = 0; i < total_points_x; i++)
	{
		for (int j = 0; j < total_points_y; j++)
		{
			if (marker[i][j] == 0) continue;

			qField[IR][i][j] = 1.0;
			qField[IU][i][j] = 3.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = 0.71429;
		}
	}
}

void GLOBAL::Compute_Boundary()
{
	vector< vector< int > >& marker = mesh->Get_Marker();
	//左边界：超声速入口
	for (int i = 0; i < num_ghost_point; i++)
	{
		for (int j = num_ghost_point; j < num_ghost_point + num_grid_point_y; j++)
		{
			qField[IR][i][j] = 1.0;
			qField[IU][i][j] = 3.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = 0.71429;
		}
	}

	//上边界：outflow
	for (int i = num_ghost_point; i < num_ghost_point + num_grid_point_x; i++)
	{
		for (int j = num_ghost_point + num_grid_point_y; j < total_points_y; j++)
		{
			qField[IR][i][j] = qField[IR][i][j - 1];
			qField[IU][i][j] = qField[IU][i][j - 1];
			qField[IV][i][j] = qField[IV][i][j - 1];
			qField[IP][i][j] = qField[IP][i][j - 1];
		}
	}

	//右边界：outflow
	for (int i = num_ghost_point + num_grid_point_x; i < total_points_x; i++)
	{
		for (int j = num_ghost_point; j < num_ghost_point + num_grid_point_y; j++)
		{
			if (marker[i - 1][j] == 0) continue;

			qField[IR][i][j] = qField[IR][i - 1][j];
			qField[IU][i][j] = qField[IU][i - 1][j];
			qField[IV][i][j] = qField[IV][i - 1][j];
			qField[IP][i][j] = qField[IP][i - 1][j];
		}
	}

	//下边界, no-slip wall
	for (int i = num_ghost_point; i < num_ghost_point + num_grid_point_x; i++)
	{
		for (int j = num_ghost_point - 1; j >= 0; j--)
		{
			qField[IR][i][j] = qField[IR][i][j + 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j + 1];
		}
	}

	//左物面, no-slip wall
	for (int i = num_ghost_point + Iw; i < num_ghost_point + Iw + num_ghost_point; i++)
	{
		for (int j = num_ghost_point + Jw1; j < num_ghost_point + Jw2; j++)
		{
			qField[IR][i][j] = qField[IR][i - 1][j];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i - 1][j];
		}
	}

	//上物面, no-slip wall
	for (int i = num_ghost_point + Iw; i < num_ghost_point + num_grid_point_x; i++)
	{
		//for (int j = Jw2; j < num_ghost_point + Jw2; j++)
		for (int j = num_ghost_point + Jw2 - 1; j >= Jw2; j--)
		{
			qField[IR][i][j] = qField[IR][i][j + 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j + 1];
		}
	}

	//下物面, no-slip wall
	for (int i = num_ghost_point + Iw; i < num_ghost_point + num_grid_point_x; i++)
	{
		for (int j = num_ghost_point + Jw1; j < num_ghost_point + Jw1 + num_ghost_point; j++)
		{
			qField[IR][i][j] = qField[IR][i][j - 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j - 1];
		}
	}
}