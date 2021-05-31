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
string tec_file_name;

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
	flow_save_steps		  = 200;	//流场输出间隔步数
	converge_criterion	  = 1e-8;	//残差收敛标准
	tec_file_name		  = "../../flow.plt";
}

void Init_Flow()
{
	//流场初始化
	Allocate_3D_Vector(qField,	  total_points_x, total_points_y, num_of_prim_vars);
	Allocate_3D_Vector(qField_N1, total_points_x, total_points_y, num_of_prim_vars);
	//Allocate_3D_Vector(qField_N2, total_points_x, total_points_y, num_of_prim_vars);
	//Allocate_3D_Vector(qField_N3, total_points_x, total_points_y, num_of_prim_vars);

	//流场赋初值
	VInt2D& marker = mesh->Get_Marker();

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			if (marker[i][j] == 0) continue;

			qField[i][j][IR] = 1.0;
			qField[i][j][IU] = 3.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = 0.71429;
		}
	}

	qField_N1 = qField;
}

void Compute_Boundary()
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

void Load_Q()
{
	qField = qField_N1;
}

void Set_Solve_Direction(char direction)
{
	solve_direction = direction;
}

void Primitive_To_Conservative( VDouble & primitive, VDouble &conservative )
{
	double gama = 1.4;
	double rho, u, v, p;
	rho = primitive[IR];
	u = primitive[IU];
	v = primitive[IV];
	p = primitive[IP];

	conservative[IR] = rho;
	conservative[IU] = rho * u;
	conservative[IV] = rho * v;
	conservative[IP] = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
}

void Conservative_To_Primitive(VDouble& conservative, VDouble &primitive)
{
	double gama = 1.4;
	double rho		= conservative[IR];
	double rho_u	= conservative[IU];
	double rho_v	= conservative[IV];
	double E		= conservative[IP];

	double u = rho_u / rho;
	double v = rho_v / rho;
	double p = (gama - 1) * (E - 0.5 * rho * (u * u + v * v));

	primitive[IR] = rho;
	primitive[IU] = u;
	primitive[IV] = v;
	primitive[IP] = p;
}

//计算点的标号范围
void Get_IJK_Region(int& ist, int& ied, int& jst, int& jed)
{
	ist = num_ghost_point;
	ied = num_ghost_point + num_grid_point_x;

	jst = num_ghost_point;
	jed = num_ghost_point + num_grid_point_y;
}