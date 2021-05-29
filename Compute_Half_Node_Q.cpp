#include <iostream>
#include "Compute_Half_Node_Q.h"
#include "Global_Variables.h"
#include "Geometry.h"
#include "2D_Euler_Solver.h"

using namespace GLOBAL;

void GLOBAL::Half_Node_Q()
{
	auto * half_node_q = new Half_Node_Q_Solver();

	half_node_q->Half_Node_Q();

	delete half_node_q;
}

Half_Node_Q_Solver::Half_Node_Q_Solver()
{
	half_node_Q_l.resize(num_of_prim_vars);
	half_node_Q_r.resize(num_of_prim_vars);

	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		Allocate_2D_Vector(half_node_Q_l[iVar], num_half_point_x, num_half_point_y);
		Allocate_2D_Vector(half_node_Q_r[iVar], num_half_point_x, num_half_point_y);
	}
}

void Half_Node_Q_Solver::Half_Node_Q()
{
	if (method_of_half_q == 1)
	{
		this->Half_Node_Q_MUSCL();
	}
	else if (method_of_half_q == 2)
	{
		//WENO无需进行半节点插值
	}
	else if (method_of_half_q == 3)
	{
		this->Half_Node_Q_WCNS();
	}
	else
	{
		cout << "半节点变量值插值方法错误，请检查！" << endl;
	}
}

void Half_Node_Q_Solver::Half_Node_Q_MUSCL()
{
	//先在x方向进行插值
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		vector< vector< double > >& qField0 = qField[iVar];
		for (int j = 0; j < num_half_point_y; j++)
		{
			for (int i = 0; i < num_half_point_x; i++)
			{
				double du_p1 = qField0[i + 1][j] - qField0[i    ][j];
				double du_m1 = qField0[i    ][j] - qField0[i - 1][j];
				double du_p3 = qField0[i + 2][j] - qField0[i + 1][j];

				double ita_m1_p = du_p1 / du_m1;
				double ita_p1_m = du_m1 / du_p1;
				double ita_p3_m = du_p1 / du_p3;
				double ita_p1_p = du_p3 / du_p1;
				
				double fai1 = Limiter_Function(ita_m1_p);
				double fai2 = Limiter_Function(ita_p1_m);
				double fai3 = Limiter_Function(ita_p3_m);
				double fai4 = Limiter_Function(ita_p1_p);

				half_node_Q_l[iVar][i][j] = qField0[i    ][j] + 1.0 / 4.0 * ((1 - muscl_k) * fai1 * du_m1 
																		   + (1 + muscl_k) * fai2 * du_p1);

				half_node_Q_r[iVar][i][j] = qField0[i + 1][j] - 1.0 / 4.0 * ((1 - muscl_k) * fai3 * du_p3 
																		   + (1 + muscl_k) * fai4 * du_p1);
			}
		}
	}
}

void Half_Node_Q_Solver::Half_Node_Q_MUSCL_Y()
{
	//在y方向进行插值
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		vector< vector< double > >& qField0 = qField[iVar];
		for (int i = 0; i < num_half_point_x; i++)
		{
			for (int j = 0; j < num_half_point_y; j++)
			{
				double du_p1 = qField0[i][j + 1] - qField0[i][j    ];
				double du_m1 = qField0[i][j    ] - qField0[i][j - 1];
				double du_p3 = qField0[i][j + 2] - qField0[i][j + 1];

				double ita_m1_p = du_p1 / du_m1;
				double ita_p1_m = du_m1 / du_p1;
				double ita_p3_m = du_p1 / du_p3;
				double ita_p1_p = du_p3 / du_p1;
				
				double fai1 = Limiter_Function(ita_m1_p);
				double fai2 = Limiter_Function(ita_p1_m);
				double fai3 = Limiter_Function(ita_p3_m);
				double fai4 = Limiter_Function(ita_p1_p);

				half_node_Q_l[iVar][i][j] = qField0[i][j	] + 1.0 / 4.0 * ((1 - muscl_k) * fai1 * du_m1 
																		   + (1 + muscl_k) * fai2 * du_p1);

				half_node_Q_r[iVar][i][j] = qField0[i][j + 1] - 1.0 / 4.0 * ((1 - muscl_k) * fai3 * du_p3 
																		   + (1 + muscl_k) * fai4 * du_p1);
			}
		}
	}
}

void Half_Node_Q_Solver::Half_Node_Q_WCNS()
{

}

double Half_Node_Q_Solver::Limiter_Function( double ita )
{
	return vanleer_limiter(ita, 1.0);
}

//限制器函数
double minmod_limiter(double a, double b)
{
	if (a * b <= 0)
		return 0;
	else
	{
		if ((fabs(a) - fabs(b)) > 0)
			return b;
		else
			return a;
	}
}

double vanleer_limiter(double a, double)
{
	return (a + abs(a)) / (1.0 + abs(a));
}

double superbee_limiter(double a, double)
{
	double tmp1 = min(2.0 * a, 1.0);
	double tmp2 = min(a, 2.0);

	return max(0.0, max(tmp1, tmp2));
}
