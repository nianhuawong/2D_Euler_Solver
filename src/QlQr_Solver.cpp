#include <iostream>
#include "QlQr_Solver.h"
#include "Geometry.h"
#include "2D_Euler_Solver.h"

void Solve_QlQr()
{
	auto * half_node_q = new QlQr_Solver();

	half_node_q->Solve_QlQr();

	delete half_node_q;
}

VDouble3D qField;
VDouble3D qField1;
VDouble3D qField2;
VDouble3D qField_N0;
VDouble3D qField_N1;

QlQr_Solver::QlQr_Solver()
{
	Get_IJK_Region(ist, ied, jst, jed);
	Allocate_3D_Vector(qField1, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(qField2, num_half_point_x, num_half_point_y, num_of_prim_vars);
}

void QlQr_Solver::Solve_QlQr()
{
	if (method_of_half_q == 1)
	{
		this->QlQr_MUSCL();
	}
	else if (method_of_half_q == 2)
	{
		//WENO无需进行半节点插值
	}
	else if (method_of_half_q == 3)
	{
		this->QlQr_WCNS();
	}
	else
	{
		cout << "半节点变量值插值方法错误，请检查！" << endl;
	}
}

void QlQr_Solver::QlQr_MUSCL()
{
	if (solve_direction == 'x')
	{
		this->QlQr_MUSCL_X();
	}
	else if (solve_direction == 'y')
	{
		this->QlQr_MUSCL_Y();
	}
	else
	{
		cout << "MUSCL插值出错，请检查！" << endl;
	}
}

void QlQr_Solver::QlQr_MUSCL_X()
{
	//在x方向进行插值
	VInt2D& marker = mesh->Get_Marker();
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;

			VDouble qVector_m1 = qField[i - 1][j];
			VDouble qVector_c0 = qField[i    ][j];
			VDouble qVector_p1 = qField[i + 1][j];
			VDouble qVector_p2 = qField[i + 2][j];

			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				double du_p1 = qVector_p1[iVar] - qVector_c0[iVar] + SMALL;
				double du_m1 = qVector_c0[iVar] - qVector_m1[iVar] + SMALL;
				double du_p3 = qVector_p2[iVar] - qVector_p1[iVar] + SMALL;

				double ita_m1_p = du_p1 / du_m1;
				double ita_p1_m = du_m1 / du_p1;
				double ita_p3_m = du_p1 / du_p3;
				double ita_p1_p = du_p3 / du_p1;

				double fai1 = Limiter_Function(ita_m1_p);
				double fai2 = Limiter_Function(ita_p1_m);
				double fai3 = Limiter_Function(ita_p3_m);
				double fai4 = Limiter_Function(ita_p1_p);

				qField1[i][j][iVar] = qVector_c0[iVar] + 1.0 / 4.0 * ((1 - muscl_k) * fai1 * du_m1
																	+ (1 + muscl_k) * fai2 * du_p1);

				qField2[i][j][iVar] = qVector_p1[iVar] - 1.0 / 4.0 * ((1 - muscl_k) * fai3 * du_p3
																	+ (1 + muscl_k) * fai4 * du_p1);
				//if (IsNaN(qField1[i][j]) || IsNaN(qField2[i][j]))
				//{
				//	int kkk = 1;
				//}
			}
		}
	}
}

void QlQr_Solver::QlQr_MUSCL_Y()
{
	//在y方向进行插值
	VInt2D& marker = mesh->Get_Marker();
	for (int i = ist; i < ied - 1; i++)
	{
		for (int j = jst; j < jed - 1; j++)
		{
			if (marker[i][j] == 0) continue;

			VDouble qVector_m1 = qField[i][j - 1];
			VDouble qVector_c0 = qField[i][j    ];
			VDouble qVector_p1 = qField[i][j + 1];
			VDouble qVector_p2 = qField[i][j + 2];

			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				double du_p1 = qVector_p1[iVar] - qVector_c0[iVar] + SMALL;
				double du_m1 = qVector_c0[iVar] - qVector_m1[iVar] + SMALL;
				double du_p3 = qVector_p2[iVar] - qVector_p1[iVar] + SMALL;

				double ita_m1_p = du_p1 / du_m1;
				double ita_p1_m = du_m1 / du_p1;
				double ita_p3_m = du_p1 / du_p3;
				double ita_p1_p = du_p3 / du_p1;

				double fai1 = Limiter_Function(ita_m1_p);
				double fai2 = Limiter_Function(ita_p1_m);
				double fai3 = Limiter_Function(ita_p3_m);
				double fai4 = Limiter_Function(ita_p1_p);

				qField1[i][j][iVar] = qVector_c0[iVar] + 1.0 / 4.0 * ((1 - muscl_k) * fai1 * du_m1
																	+ (1 + muscl_k) * fai2 * du_p1);

				qField2[i][j][iVar] = qVector_p1[iVar] - 1.0 / 4.0 * ((1 - muscl_k) * fai3 * du_p3
																	+ (1 + muscl_k) * fai4 * du_p1);
			}
		}
	}
}

void QlQr_Solver::QlQr_WCNS()
{
	if (solve_direction == 'x')
	{
		this->QlQr_WCNS_X();
	}
	else if (solve_direction == 'y')
	{
		this->QlQr_WCNS_Y();
	}
}

void QlQr_Solver::QlQr_WCNS_X()
{
	//计算Lagrange插值系数
	VDouble3D g1, g2, g3;
	Allocate_3D_Vector(g1, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(g2, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(g3, num_half_point_x, num_half_point_y, num_of_prim_vars);

	VDouble3D s1, s2, s3;
	Allocate_3D_Vector(s1, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(s2, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(s3, num_half_point_x, num_half_point_y, num_of_prim_vars);

	double ds = dx;

	VInt2D& marker = mesh->Get_Marker();
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				g1[i][j][iVar] = (		 qField[i - 2][j][iVar] - 4.0 * qField[i - 1][j][iVar] + 3.0 *	qField[i    ][j][iVar]) / 2.0 / ds;
				g2[i][j][iVar] = (		 qField[i + 1][j][iVar] -		qField[i - 1][j][iVar]								  ) / 2.0 / ds;
				g3[i][j][iVar] = (-3.0 * qField[i    ][j][iVar] + 4.0 * qField[i + 1][j][iVar] -		qField[i + 2][j][iVar]) / 2.0 / ds;

				s1[i][j][iVar] = (		 qField[i - 2][j][iVar] - 2.0 * qField[i - 1][j][iVar] +		qField[i    ][j][iVar]) / ds / ds;
				s2[i][j][iVar] = (		 qField[i - 1][j][iVar] - 2.0 * qField[i    ][j][iVar] +		qField[i + 1][j][iVar]) / ds / ds;
				s3[i][j][iVar] = (		 qField[i    ][j][iVar] - 2.0 * qField[i + 1][j][iVar] +		qField[i + 2][j][iVar]) / ds / ds;
			}
		}
	}
	//计算光滑因子
	VDouble3D IS1, IS2, IS3;
	Allocate_3D_Vector(IS1, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS2, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS3, num_half_point_x, num_half_point_y, num_of_prim_vars);

	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				IS1[i][j][iVar] = pow((ds * g1[i][j][iVar]), 2) + pow((ds * ds * s1[i][j][iVar]), 2);
				IS2[i][j][iVar] = pow((ds * g2[i][j][iVar]), 2) + pow((ds * ds * s2[i][j][iVar]), 2);
				IS3[i][j][iVar] = pow((ds * g3[i][j][iVar]), 2) + pow((ds * ds * s3[i][j][iVar]), 2);
			}
		}
	}

	double C11 = 1.0 / 16.0, C21 = 10.0 / 16.0, C31 = 5.0 / 16.0;
	double C12 = 5.0 / 16.0, C22 = 10.0 / 16.0, C32 = 1.0 / 16.0;
	double eps = 1e-6;

	//j+1/2处的变量左右值和通量
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				//j+1/2处的左右值，可以由3个模板插值得到3个值
				double qField_Left1 = qField[i][j][iVar] + ds * g1[i][j][iVar] / 2.0 + ds * ds * s1[i][j][iVar] / 8.0;
				double qField_Left2 = qField[i][j][iVar] + ds * g2[i][j][iVar] / 2.0 + ds * ds * s2[i][j][iVar] / 8.0;
				double qField_Left3 = qField[i][j][iVar] + ds * g3[i][j][iVar] / 2.0 + ds * ds * s3[i][j][iVar] / 8.0;

				double qField_Right1 = qField[i + 1][j][iVar] - ds * g1[i + 1][j][iVar] / 2.0 + ds * ds * s1[i + 1][j][iVar] / 8.0;
				double qField_Right2 = qField[i + 1][j][iVar] - ds * g2[i + 1][j][iVar] / 2.0 + ds * ds * s2[i + 1][j][iVar] / 8.0;
				double qField_Right3 = qField[i + 1][j][iVar] - ds * g3[i + 1][j][iVar] / 2.0 + ds * ds * s3[i + 1][j][iVar] / 8.0;

				//取非线性加权，左右值分别取不同的权值
				double a1 = C11 / pow((eps + IS1[i][j][iVar]), 2);
				double a2 = C21 / pow((eps + IS2[i][j][iVar]), 2);
				double a3 = C31 / pow((eps + IS3[i][j][iVar]), 2);

				double w1 = a1 / (a1 + a2 + a3);
				double w2 = a2 / (a1 + a2 + a3);
				double w3 = a3 / (a1 + a2 + a3);

				qField1[i][j][iVar] = w1 * qField_Left1 + w2 * qField_Left2 + w3 * qField_Left3;

				//取非线性加权，左右值分别取不同的权值
				a1 = C12 / pow((eps + IS1[i][j][iVar]), 2);
				a2 = C22 / pow((eps + IS2[i][j][iVar]), 2);
				a3 = C32 / pow((eps + IS3[i][j][iVar]), 2);

				w1 = a1 / (a1 + a2 + a3);
				w2 = a2 / (a1 + a2 + a3);
				w3 = a3 / (a1 + a2 + a3);
				qField2[i][j][iVar] = w1 * qField_Right1 + w2 * qField_Right2 + w3 * qField_Right3;
			}
		}
	}
}

void QlQr_Solver::QlQr_WCNS_Y()
{
	//计算Lagrange插值系数
	VDouble3D g1, g2, g3;
	Allocate_3D_Vector(g1, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(g2, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(g3, num_half_point_x, num_half_point_y, num_of_prim_vars);

	VDouble3D s1, s2, s3;
	Allocate_3D_Vector(s1, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(s2, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(s3, num_half_point_x, num_half_point_y, num_of_prim_vars);

	double ds = dy;

	VInt2D& marker = mesh->Get_Marker();
	for (int i = ist; i < ied - 1; i++)
	{
		for (int j = jst; j < jed - 1; j++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				g1[i][j][iVar] = (		 qField[i][j - 2][iVar] - 4.0 * qField[i][j - 1][iVar] + 3.0 *	qField[i][j    ][iVar]) / 2.0 / ds;
				g2[i][j][iVar] = (		 qField[i][j + 1][iVar] -		qField[i][j - 1][iVar]								  ) / 2.0 / ds;
				g3[i][j][iVar] = (-3.0 * qField[i][j    ][iVar] + 4.0 * qField[i][j + 1][iVar] -		qField[i][j + 2][iVar]) / 2.0 / ds;

				s1[i][j][iVar] = (		 qField[i][j - 2][iVar] - 2.0 * qField[i][j - 1][iVar] +		qField[i][j    ][iVar]) / ds / ds;
				s2[i][j][iVar] = (		 qField[i][j - 1][iVar] - 2.0 * qField[i][j    ][iVar] +		qField[i][j + 1][iVar]) / ds / ds;
				s3[i][j][iVar] = (		 qField[i][j    ][iVar] - 2.0 * qField[i][j + 1][iVar] +		qField[i][j + 2][iVar]) / ds / ds;
			}
		}
	}
	//计算光滑因子
	VDouble3D IS1, IS2, IS3;
	Allocate_3D_Vector(IS1, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS2, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS3, num_half_point_x, num_half_point_y, num_of_prim_vars);

	for (int i = ist; i < ied - 1; i++)
	{
		for (int j = jst; j < jed - 1; j++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				IS1[i][j][iVar] = pow((ds * g1[i][j][iVar]), 2) + pow((ds * ds * s1[i][j][iVar]), 2);
				IS2[i][j][iVar] = pow((ds * g2[i][j][iVar]), 2) + pow((ds * ds * s2[i][j][iVar]), 2);
				IS3[i][j][iVar] = pow((ds * g3[i][j][iVar]), 2) + pow((ds * ds * s3[i][j][iVar]), 2);
			}
		}
	}

	double C11 = 1.0 / 16.0, C21 = 10.0 / 16.0, C31 = 5.0 / 16.0;
	double C12 = 5.0 / 16.0, C22 = 10.0 / 16.0, C32 = 1.0 / 16.0;
	double eps = 1e-6;

	//j+1/2处的变量左右值和通量
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				//j+1/2处的左右值，可以由3个模板插值得到3个值
				double qField_Left1 = qField[i][j][iVar] + ds * g1[i][j][iVar] / 2.0 + ds * ds * s1[i][j][iVar] / 8.0;
				double qField_Left2 = qField[i][j][iVar] + ds * g2[i][j][iVar] / 2.0 + ds * ds * s2[i][j][iVar] / 8.0;
				double qField_Left3 = qField[i][j][iVar] + ds * g3[i][j][iVar] / 2.0 + ds * ds * s3[i][j][iVar] / 8.0;

				double qField_Right1 = qField[i][j + 1][iVar] - ds * g1[i][j + 1][iVar] / 2.0 + ds * ds * s1[i][j + 1][iVar] / 8.0;
				double qField_Right2 = qField[i][j + 1][iVar] - ds * g2[i][j + 1][iVar] / 2.0 + ds * ds * s2[i][j + 1][iVar] / 8.0;
				double qField_Right3 = qField[i][j + 1][iVar] - ds * g3[i][j + 1][iVar] / 2.0 + ds * ds * s3[i][j + 1][iVar] / 8.0;

				//取非线性加权，左右值分别取不同的权值
				double a1 = C11 / pow((eps + IS1[i][j][iVar]), 2);
				double a2 = C21 / pow((eps + IS2[i][j][iVar]), 2);
				double a3 = C31 / pow((eps + IS3[i][j][iVar]), 2);

				double w1 = a1 / (a1 + a2 + a3);
				double w2 = a2 / (a1 + a2 + a3);
				double w3 = a3 / (a1 + a2 + a3);

				qField1[i][j][iVar] = w1 * qField_Left1 + w2 * qField_Left2 + w3 * qField_Left3;

				//取非线性加权，左右值分别取不同的权值
				a1 = C12 / pow((eps + IS1[i][j][iVar]), 2);
				a2 = C22 / pow((eps + IS2[i][j][iVar]), 2);
				a3 = C32 / pow((eps + IS3[i][j][iVar]), 2);

				w1 = a1 / (a1 + a2 + a3);
				w2 = a2 / (a1 + a2 + a3);
				w3 = a3 / (a1 + a2 + a3);
				qField2[i][j][iVar] = w1 * qField_Right1 + w2 * qField_Right2 + w3 * qField_Right3;
			}
		}
	}
}


double QlQr_Solver::Limiter_Function( double ita )
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
