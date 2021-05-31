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
VDouble3D qField_N1;
VDouble3D qField_N2;
VDouble3D qField_N3;
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
	vector< vector< int > >& marker = mesh->Get_Marker();
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
			}
		}
	}
}

void QlQr_Solver::QlQr_MUSCL_Y()
{
	//在y方向进行插值
	vector< vector< int > >& marker = mesh->Get_Marker();
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
