#include <iostream>
#include <cmath>
#include "Global.h"
#include "QlQr_Solver.h"
#include "Flux_Solver.h"
#include "2D_Euler_Solver.h"
#include "Geometry.h"

VDouble3D fluxVector;

void Solve_Flux()
{
	Flux_Solver* fluxSolver = new Flux_Solver();

	fluxSolver->Solve_Flux();

	delete fluxSolver;
}

Flux_Solver::Flux_Solver()
{
	this->gama = 1.4;

	Get_IJK_Region(ist, ied, jst, jed);

	Allocate_3D_Vector(fluxVector,  num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(fluxVector1, total_points_x, total_points_y, num_of_prim_vars);
	Allocate_3D_Vector(fluxVector2, total_points_x, total_points_y, num_of_prim_vars);

}

void Flux_Solver::Solve_Flux()
{
	if (method_of_flux == 1)
	{
		Roe_Scheme();
	}
	else if (method_of_flux == 2)
	{	
		Steger_Warming_Scheme();//二阶Steger_Warming
	}
	else if (method_of_flux == 3)
	{
		VanLeer_Scheme();		//二阶vanLeer
	}
	else if (method_of_flux == 4)
	{
		//半节点左右通量插值：WENO方法
		WENO_Scheme();		
	}
	else if (method_of_flux == 5)
	{
		//WCNS原始变量插值已经做完；
		//计算半节点左右通量和roe通量
		Roe_Scheme();	
	}
	else
	{
		cout << "flux计算方法选择错误，请检查！" << endl;
	}
}

void Flux_Solver::WENO_Scheme()
{
	if (solve_direction == 'x')
	{
		this->Steger_Warming_Scheme_X();  //整数节点通量计算：Steger-Warming方法
		this->WENO_Scheme_X();
	}
	else if (solve_direction == 'y')
	{
		this->Steger_Warming_Scheme_Y();
		this->WENO_Scheme_Y();
	}
	else
	{
		cout << "出错，请检查！" << endl;
	}
}

void Flux_Solver::WENO_Scheme_X()
{
	//q(j+1/2)-
	VDouble3D q11, q21, q31;
	Allocate_3D_Vector(q11, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(q21, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(q31, num_half_point_x, num_half_point_y, num_of_prim_vars);

	VInt2D& marker = mesh->Get_Marker();
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 1; i <= ied - 1 ; i++)//i=0, ied, ied + 1三个点没有值
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				q11[i][j][iVar] = -1.0 / 6.0 * fluxVector1[i - 1][j][iVar] + 5.0 / 6.0 * fluxVector1[i    ][j][iVar] + 1.0 / 3.0 * fluxVector1[i + 1][j][iVar];
				q21[i][j][iVar] =  1.0 / 3.0 * fluxVector1[i    ][j][iVar] + 5.0 / 6.0 * fluxVector1[i + 1][j][iVar] - 1.0 / 6.0 * fluxVector1[i + 2][j][iVar];
				q31[i][j][iVar] = 11.0 / 6.0 * fluxVector1[i + 1][j][iVar] - 7.0 / 6.0 * fluxVector1[i + 2][j][iVar] + 1.0 / 3.0 * fluxVector1[i + 3][j][iVar];
			}
		}
	}

	//q(j+1/2)+
	VDouble3D q12, q22, q32;
	Allocate_3D_Vector(q12, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(q22, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(q32, num_half_point_x, num_half_point_y, num_of_prim_vars);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 2; i <= ied; i++)//i=0, 1, ied + 1三个点没有值
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				q12[i][j][iVar] =  1.0 / 3.0 * fluxVector2[i - 2][j][iVar] - 7.0 / 6.0 * fluxVector2[i - 1][j][iVar] + 11.0 / 6.0 * fluxVector2[i    ][j][iVar];
				q22[i][j][iVar] = -1.0 / 6.0 * fluxVector2[i - 1][j][iVar] + 5.0 / 6.0 * fluxVector2[i    ][j][iVar] +  1.0 / 3.0 * fluxVector2[i + 1][j][iVar];
				q32[i][j][iVar] =  1.0 / 3.0 * fluxVector2[i    ][j][iVar] + 5.0 / 6.0 * fluxVector2[i + 1][j][iVar] -  1.0 / 6.0 * fluxVector2[i + 2][j][iVar];
			}
		}
	}

	//IS(j+1/2)-
	VDouble3D IS11, IS21, IS31;
	Allocate_3D_Vector(IS11, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS21, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS31, num_half_point_x, num_half_point_y, num_of_prim_vars);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 1; i <= ied - 1; i++)//i=0, ied, ied + 1三个点没有值
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				IS11[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector1[i - 1][j][iVar] - 2.0 * fluxVector1[i    ][j][iVar] +       fluxVector1[i + 1][j][iVar]), 2)
								   + 1.0 / 4.0 * pow((fluxVector1[i - 1][j][iVar] - 4.0 * fluxVector1[i    ][j][iVar] + 3.0 * fluxVector1[i + 1][j][iVar]), 2);

				IS21[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector1[i    ][j][iVar] - 2.0 * fluxVector1[i + 1][j][iVar] +       fluxVector1[i + 2][j][iVar]), 2)
								   + 1.0 / 4.0 * pow((fluxVector1[i    ][j][iVar] -       fluxVector1[i + 2][j][iVar])									  , 2);

				IS31[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector1[i + 1][j][iVar] - 2.0 * fluxVector1[i + 2][j][iVar] +		  fluxVector1[i + 3][j][iVar]), 2)
								   + 1.0 / 4.0 * pow((fluxVector1[i + 3][j][iVar] - 4.0 * fluxVector1[i + 2][j][iVar] + 3.0 * fluxVector1[i + 1][j][iVar]), 2);

			}
		}
	}

	//IS(j+1/2)+
	VDouble3D IS12, IS22, IS32;
	Allocate_3D_Vector(IS12, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS22, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS32, num_half_point_x, num_half_point_y, num_of_prim_vars);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 2; i <= ied; i++)//i=0, 1, ied + 1三个点没有值
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				IS12[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector2[i - 2][j][iVar] - 2.0 * fluxVector2[i - 1][j][iVar] +       fluxVector2[i    ][j][iVar]), 2)
								 +  1.0 / 4.0  * pow((fluxVector2[i - 2][j][iVar] - 4.0 * fluxVector2[i - 1][j][iVar] + 3.0 * fluxVector2[i    ][j][iVar]), 2);

				IS22[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector2[i - 1][j][iVar] - 2.0 * fluxVector2[i    ][j][iVar] +       fluxVector2[i + 1][j][iVar]), 2)
								 +  1.0 / 4.0  * pow((fluxVector2[i - 1][j][iVar] -       fluxVector2[i + 1][j][iVar])									  , 2);

				IS32[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector2[i    ][j][iVar] - 2.0 * fluxVector2[i + 1][j][iVar] +       fluxVector2[i + 2][j][iVar]), 2)
								 +  1.0 / 4.0  * pow((fluxVector2[i + 2][j][iVar] - 4.0 * fluxVector2[i + 1][j][iVar] + 3.0 * fluxVector2[i    ][j][iVar]), 2);
			}
		}
	}

	double C11 = 3.0 / 10.0, C21 = 3.0 / 5.0, C31 = 1.0 / 10.0;
	double C12 = 1.0 / 10.0, C22 = 3.0 / 5.0, C32 = 3.0 / 10.0;
	double eps = 1e-6;

	//j+1/2(-)处的通量
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 2; i <= ied - 1; i++)//i = 0, 1, ied, ied + 1没有计算，在后续计算
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				double a1 = C11 / pow((eps + IS11[i][j][iVar]), 2);
				double a2 = C21 / pow((eps + IS21[i][j][iVar]), 2);
				double a3 = C31 / pow((eps + IS31[i][j][iVar]), 2);

				double w1 = a1 / (a1 + a2 + a3);
				double w2 = a2 / (a1 + a2 + a3);
				double w3 = a3 / (a1 + a2 + a3);

				fluxVector1[i][j][iVar] = w1 * q11[i][j][iVar] + w2 * q21[i][j][iVar] + w3 * q31[i][j][iVar];
			}
		}
	}

	//j+1/2(+)处的通量
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 2; i <= ied - 1; i++)//i = 0, 1, ied, ied + 1没有计算，在后续计算
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				double a1 = C12 / pow((eps + IS12[i][j][iVar]), 2);
				double a2 = C22 / pow((eps + IS22[i][j][iVar]), 2);
				double a3 = C32 / pow((eps + IS32[i][j][iVar]), 2);

				double w1 = a1 / (a1 + a2 + a3);
				double w2 = a2 / (a1 + a2 + a3);
				double w3 = a3 / (a1 + a2 + a3);

				fluxVector2[i][j][iVar] = w1 * q12[i][j][iVar] + w2 * q22[i][j][iVar] + w3 * q32[i][j][iVar];
			}
		}
	}

	//i=0,1,ied,ied+1处的通量值
	for (int j = jst; j <= jed; j++)
	{
		for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
		{
			fluxVector1[0][j][iVar] = fluxVector1[2][j][iVar];
			fluxVector1[1][j][iVar] = fluxVector1[2][j][iVar];
			fluxVector1[ied    ][j][iVar] = fluxVector1[ied - 1][j][iVar];
			fluxVector1[ied + 1][j][iVar] = fluxVector1[ied - 1][j][iVar];

			fluxVector2[0][j][iVar] = fluxVector2[2][j][iVar];
			fluxVector2[1][j][iVar] = fluxVector2[2][j][iVar];
			fluxVector2[ied    ][j][iVar] = fluxVector2[ied - 1][j][iVar];
			fluxVector2[ied + 1][j][iVar] = fluxVector2[ied - 1][j][iVar];
		}
	}

	//半节点通量值：Steger_Warming进行正负通量相加
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 0; i <= ied + 1; i++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				fluxVector[i][j][iVar] = fluxVector1[i][j][iVar] + fluxVector2[i][j][iVar];
			}
		}
	}
}

void Flux_Solver::WENO_Scheme_Y()
{
	//q(j+1/2)-
	VDouble3D q11, q21, q31;
	Allocate_3D_Vector(q11, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(q21, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(q31, num_half_point_x, num_half_point_y, num_of_prim_vars);

	VInt2D& marker = mesh->Get_Marker();
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 1; j <= jed - 1; j++)//j=0, jed, jed + 1三个点没有值
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				q11[i][j][iVar] = -1.0 / 6.0 * fluxVector1[i][j - 1][iVar] + 5.0 / 6.0 * fluxVector1[i][j    ][iVar] + 1.0 / 3.0 * fluxVector1[i][j + 1][iVar];
				q21[i][j][iVar] =  1.0 / 3.0 * fluxVector1[i][j    ][iVar] + 5.0 / 6.0 * fluxVector1[i][j + 1][iVar] - 1.0 / 6.0 * fluxVector1[i][j + 2][iVar];
				q31[i][j][iVar] = 11.0 / 6.0 * fluxVector1[i][j + 1][iVar] - 7.0 / 6.0 * fluxVector1[i][j + 2][iVar] + 1.0 / 3.0 * fluxVector1[i][j + 3][iVar];
			}
		}
	}

	//q(j+1/2)+
	VDouble3D q12, q22, q32;
	Allocate_3D_Vector(q12, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(q22, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(q32, num_half_point_x, num_half_point_y, num_of_prim_vars);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 2; j <= jed; j++)//i=0, 1, ied + 1三个点没有值
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				q12[i][j][iVar] =  1.0 / 3.0 * fluxVector2[i][j - 2][iVar] - 7.0 / 6.0 * fluxVector2[i][j - 1][iVar] + 11.0 / 6.0 * fluxVector2[i][j    ][iVar];
				q22[i][j][iVar] = -1.0 / 6.0 * fluxVector2[i][j - 1][iVar] + 5.0 / 6.0 * fluxVector2[i][j    ][iVar] +  1.0 / 3.0 * fluxVector2[i][j + 1][iVar];
				q32[i][j][iVar] =  1.0 / 3.0 * fluxVector2[i][j    ][iVar] + 5.0 / 6.0 * fluxVector2[i][j + 1][iVar] -  1.0 / 6.0 * fluxVector2[i][j + 2][iVar];
			}
		}
	}

	//IS(j+1/2)-
	VDouble3D IS11, IS21, IS31;
	Allocate_3D_Vector(IS11, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS21, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS31, num_half_point_x, num_half_point_y, num_of_prim_vars);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 1; j <= jed - 1; j++)//j=0, jed, jed + 1三个点没有值
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				IS11[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector1[i][j - 1][iVar] - 2.0 * fluxVector1[i][j    ][iVar] +       fluxVector1[i][j + 1][iVar]), 2)
								   + 1.0 / 4.0 * pow((fluxVector1[i][j - 1][iVar] - 4.0 * fluxVector1[i][j    ][iVar] + 3.0 * fluxVector1[i][j + 1][iVar]), 2);

				IS21[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector1[i][j    ][iVar] - 2.0 * fluxVector1[i][j + 1][iVar] +       fluxVector1[i][j + 2][iVar]), 2)
								   + 1.0 / 4.0 * pow((fluxVector1[i][j    ][iVar] -       fluxVector1[i][j + 2][iVar])									  , 2);

				IS31[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector1[i][j + 1][iVar] - 2.0 * fluxVector1[i][j + 2][iVar] +		  fluxVector1[i][j + 3][iVar]), 2)
								   + 1.0 / 4.0 * pow((fluxVector1[i][j + 3][iVar] - 4.0 * fluxVector1[i][j + 2][iVar] + 3.0 * fluxVector1[i][j + 1][iVar]), 2);

			}
		}
	}

	//IS(j+1/2)+
	VDouble3D IS12, IS22, IS32;
	Allocate_3D_Vector(IS12, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS22, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(IS32, num_half_point_x, num_half_point_y, num_of_prim_vars);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 2; j <= jed; j++)//i=0, 1, ied + 1三个点没有值
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				IS12[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector2[i][j - 2][iVar] - 2.0 * fluxVector2[i][j - 1][iVar] +       fluxVector2[i][j    ][iVar]), 2)
								 +  1.0 / 4.0  * pow((fluxVector2[i][j - 2][iVar] - 4.0 * fluxVector2[i][j - 1][iVar] + 3.0 * fluxVector2[i][j    ][iVar]), 2);

				IS22[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector2[i][j - 1][iVar] - 2.0 * fluxVector2[i][j    ][iVar] +       fluxVector2[i][j + 1][iVar]), 2)
								 +  1.0 / 4.0  * pow((fluxVector2[i][j - 1][iVar] -       fluxVector2[i][j + 1][iVar])									  , 2);

				IS32[i][j][iVar] = 13.0 / 12.0 * pow((fluxVector2[i][j    ][iVar] - 2.0 * fluxVector2[i][j + 1][iVar] +       fluxVector2[i][j + 2][iVar]), 2)
								 +  1.0 / 4.0  * pow((fluxVector2[i][j + 2][iVar] - 4.0 * fluxVector2[i][j + 1][iVar] + 3.0 * fluxVector2[i][j    ][iVar]), 2);
			}
		}
	}

	double C11 = 3.0 / 10.0, C21 = 3.0 / 5.0, C31 = 1.0 / 10.0;
	double C12 = 1.0 / 10.0, C22 = 3.0 / 5.0, C32 = 3.0 / 10.0;
	double eps = 1e-6;

	//j+1/2(-)处的通量
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 2; j <= jed - 1; j++)//j = 0, 1, jed, jed + 1没有计算，在后续计算
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				double a1 = C11 / pow((eps + IS11[i][j][iVar]), 2);
				double a2 = C21 / pow((eps + IS21[i][j][iVar]), 2);
				double a3 = C31 / pow((eps + IS31[i][j][iVar]), 2);

				double w1 = a1 / (a1 + a2 + a3);
				double w2 = a2 / (a1 + a2 + a3);
				double w3 = a3 / (a1 + a2 + a3);

				fluxVector1[i][j][iVar] = w1 * q11[i][j][iVar] + w2 * q21[i][j][iVar] + w3 * q31[i][j][iVar];
			}
		}
	}

	//j+1/2(+)处的通量
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 2; j <= jed - 1; j++)//j = 0, 1, jed, jed + 1没有计算，在后续计算
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				double a1 = C12 / pow((eps + IS12[i][j][iVar]), 2);
				double a2 = C22 / pow((eps + IS22[i][j][iVar]), 2);
				double a3 = C32 / pow((eps + IS32[i][j][iVar]), 2);

				double w1 = a1 / (a1 + a2 + a3);
				double w2 = a2 / (a1 + a2 + a3);
				double w3 = a3 / (a1 + a2 + a3);

				fluxVector2[i][j][iVar] = w1 * q12[i][j][iVar] + w2 * q22[i][j][iVar] + w3 * q32[i][j][iVar];
			}
		}
	}

	//j=0,1,jed,jed+1处的通量值
	for (int i = ist; i <= ied; i++)
	{
		for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
		{
			fluxVector1[i][0      ][iVar] = fluxVector1[i][2	  ][iVar];
			fluxVector1[i][1      ][iVar] = fluxVector1[i][2	  ][iVar];
			fluxVector1[i][jed    ][iVar] = fluxVector1[i][jed - 1][iVar];
			fluxVector1[i][jed + 1][iVar] = fluxVector1[i][jed - 1][iVar];

			fluxVector2[i][0	  ][iVar] = fluxVector2[i][2	  ][iVar];
			fluxVector2[i][1	  ][iVar] = fluxVector2[i][2	  ][iVar];
			fluxVector2[i][jed    ][iVar] = fluxVector2[i][jed - 1][iVar];
			fluxVector2[i][jed + 1][iVar] = fluxVector2[i][jed - 1][iVar];
		}
	}

	//半节点通量值：Steger_Warming进行正负通量相加
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 0; j <= jed + 1; j++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				fluxVector[i][j][iVar] = fluxVector1[i][j][iVar] + fluxVector2[i][j][iVar];
			}
		}
	}
}

void Flux_Solver::Roe_Scheme()
{
	if (solve_direction == 'x')
	{
		this->Roe_Scheme_X();
	}
	else if (solve_direction == 'y')
	{
		this->Roe_Scheme_Y();
	}
	else
	{
		cout << "出错，请检查！" << endl;
	}
}

void Flux_Solver::Roe_Scheme_X()
{
	VInt2D& marker = mesh->Get_Marker();

#ifndef _OPENMP
	VDouble fluxVector1(num_of_prim_vars);
	VDouble fluxVector2(num_of_prim_vars);
	VDouble2D Jacobian_A;
	Allocate_2D_Vector(Jacobian_A, num_of_prim_vars, num_of_prim_vars);
#endif

#ifdef _OPENMP
#pragma omp parallel
	{
		VDouble fluxVector1(num_of_prim_vars);
		VDouble fluxVector2(num_of_prim_vars);
		VDouble2D Jacobian_A;
		Allocate_2D_Vector(Jacobian_A, num_of_prim_vars, num_of_prim_vars);
#pragma omp for 	
#endif
	for (int i = 0; i <= ied + 1; i++)
	{
		for (int j = jst; j <= jed; j++)
		{
			if (marker[i][j] == 0) continue;

			double rho1, u1, v1, p1;
			ExtractValue(qField1[i][j], rho1, u1, v1, p1);
			double H1   = Enthalpy(rho1, u1, v1, p1, gama);

			Inviscid_Flux_F(fluxVector1, rho1, u1, v1, p1);

			double rho2, u2, v2, p2;
			ExtractValue(qField2[i][j], rho2, u2, v2, p2);
			double H2   = Enthalpy(rho2, u2, v2, p2, gama);

			Inviscid_Flux_F(fluxVector2, rho2, u2, v2, p2);

			double D = sqrt(rho2 / rho1);

			double rho_roe = sqrt(rho1 * rho2);
			double u_roe = (u1 + u2 * D) / (1 + D);
			double v_roe = (v1 + v2 * D) / (1 + D);
			double H_roe = (H1 + H2 * D) / (1 + D);
			double c2_roe = (gama - 1) * (H_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe));
			double c_roe = sqrt(fabs(c2_roe));//声速取绝对值

			Compute_Jacobian(Jacobian_A, u_roe, v_roe, c_roe, H_roe);

			VDouble qPrimitive1 = qField1[i][j];
			VDouble qPrimitive2 = qField2[i][j];

			VDouble qConservative1(num_of_prim_vars);
			VDouble qConservative2(num_of_prim_vars);
			Primitive_To_Conservative(qPrimitive1, qConservative1);
			Primitive_To_Conservative(qPrimitive2, qConservative2);

			VDouble dq(num_of_prim_vars);
			dq[IR] = qConservative2[IR] - qConservative1[IR];
			dq[IU] = qConservative2[IU] - qConservative1[IU];
			dq[IV] = qConservative2[IV] - qConservative1[IV];
			dq[IP] = qConservative2[IP] - qConservative1[IP];

			VDouble flux_add(num_of_prim_vars);
			MatrixMultiply(Jacobian_A, dq, flux_add, 4, 4);

			fluxVector[i][j][IR] = 0.5 * (fluxVector1[IR] + fluxVector2[IR] - flux_add[IR]);
			fluxVector[i][j][IU] = 0.5 * (fluxVector1[IU] + fluxVector2[IU] - flux_add[IU]);
			fluxVector[i][j][IV] = 0.5 * (fluxVector1[IV] + fluxVector2[IV] - flux_add[IV]);
			fluxVector[i][j][IP] = 0.5 * (fluxVector1[IP] + fluxVector2[IP] - flux_add[IP]);
		}
	}
#ifdef _OPENMP
}
#endif // _OPENMP
}

void Flux_Solver::Roe_Scheme_Y()
{
	VInt2D& marker = mesh->Get_Marker();
#ifndef _OPENMP
	VDouble2D Jacobian_A;
	Allocate_2D_Vector(Jacobian_A, num_of_prim_vars, num_of_prim_vars);
	VDouble fluxVector1(num_of_prim_vars);
	VDouble fluxVector2(num_of_prim_vars);
#endif

#ifdef _OPENMP
#pragma omp parallel
	{
		VDouble2D Jacobian_A;
		Allocate_2D_Vector(Jacobian_A, num_of_prim_vars, num_of_prim_vars);
		VDouble fluxVector1(num_of_prim_vars);
		VDouble fluxVector2(num_of_prim_vars);
#pragma omp for 	
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 0; j <= jed + 1; j++)
		{
			if (marker[i][j] == 0) continue;

			double rho1, u1, v1, p1;
			ExtractValue(qField1[i][j], rho1, u1, v1, p1);
			double H1   = Enthalpy(rho1, u1, v1, p1, gama);

			Inviscid_Flux_G(fluxVector1, rho1, u1, v1, p1);

			double rho2, u2, v2, p2;
			ExtractValue(qField2[i][j], rho2, u2, v2, p2);
			double H2   = Enthalpy(rho2, u2, v2, p2, gama);

			Inviscid_Flux_G(fluxVector2, rho2, u2, v2, p2);

			double D = sqrt(rho2 / rho1);

			double rho_roe = sqrt(rho1 * rho2);
			double u_roe = (u1 + u2 * D) / (1 + D);
			double v_roe = (v1 + v2 * D) / (1 + D);
			double H_roe = (H1 + H2 * D) / (1 + D);
			double c2_roe = (gama - 1) * (H_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe));
			double c_roe = sqrt(fabs(c2_roe));//声速取绝对值			

			Compute_Jacobian(Jacobian_A, u_roe, v_roe, c_roe, H_roe);

			VDouble qPrimitive1 = qField1[i][j];
			VDouble qPrimitive2 = qField2[i][j];

			VDouble qConservative1(num_of_prim_vars);
			VDouble qConservative2(num_of_prim_vars);
			Primitive_To_Conservative(qPrimitive1, qConservative1);
			Primitive_To_Conservative(qPrimitive2, qConservative2);

			VDouble dq(num_of_prim_vars);
			dq[IR] = qConservative2[IR] - qConservative1[IR];
			dq[IU] = qConservative2[IU] - qConservative1[IU];
			dq[IV] = qConservative2[IV] - qConservative1[IV];
			dq[IP] = qConservative2[IP] - qConservative1[IP];

			VDouble flux_add(num_of_prim_vars);
			MatrixMultiply(Jacobian_A, dq, flux_add, 4, 4);

			fluxVector[i][j][IR] = 0.5 * (fluxVector1[IR] + fluxVector2[IR] - flux_add[IR]);
			fluxVector[i][j][IU] = 0.5 * (fluxVector1[IU] + fluxVector2[IU] - flux_add[IU]);
			fluxVector[i][j][IV] = 0.5 * (fluxVector1[IV] + fluxVector2[IV] - flux_add[IV]);
			fluxVector[i][j][IP] = 0.5 * (fluxVector1[IP] + fluxVector2[IP] - flux_add[IP]);
		}
	}
#ifdef _OPENMP
}
#endif // _OPENMP
}

void Flux_Solver::VanLeer_Scheme()
{
	if (solve_direction == 'x')
	{
		this->VanLeer_Scheme_X();
	}
	else if (solve_direction == 'y')
	{
		this->VanLeer_Scheme_Y();
	}
}

void Flux_Solver::VanLeer_Scheme_X()
{
	double nx = 1.0, ny = 0.0;

	VInt2D& marker = mesh->Get_Marker();
#ifndef _OPENMP
	VDouble fluxVector11(num_of_prim_vars);
	VDouble fluxVector22(num_of_prim_vars);
	double rho1, u1, v1, p1;
	double rho2, u2, v2, p2;
#endif
#ifdef _OPENMP
#pragma omp parallel 
	{
		VDouble fluxVector11(num_of_prim_vars);
		VDouble fluxVector22(num_of_prim_vars);
		double rho1, u1, v1, p1;
		double rho2, u2, v2, p2;
#pragma omp	for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 0; i <= ied + 1; i++)//每个通量点的值都有了
		{
			if (marker[i][j] == 0) continue;

			ExtractValue(qField1[i][j], rho1, u1, v1, p1);

			double c1   = sqrt(fabs(gama * p1 / rho1));//声速取绝对值
			double vn1  = nx * u1 + ny * v1;
			double m1   = vn1 / c1;
			double h1   = Enthalpy(rho1, u1, v1, p1, gama);

			ExtractValue(qField2[i][j], rho2, u2, v2, p2);

			double c2   = sqrt(fabs(gama * p2 / rho2));//声速取绝对值
			double vn2  = nx * u2 + ny * v2;
			double m2   = vn2 / c2;
			double h2 = Enthalpy(rho2, u2, v2, p2, gama);

			if (m1 > 1.0)
			{
				Inviscid_Flux_F(fluxVector11, rho1, u1, v1, p1);
			}
			else if (m1 < -1.0)
			{
				fluxVector11[IR] = 0;
				fluxVector11[IU] = 0;
				fluxVector11[IV] = 0;
				fluxVector11[IP] = 0;
			}
			else
			{
				double fmass = rho1 * c1 * 0.25 * (m1 + 1) * (m1 + 1);
				double term0 = (-vn1 + 2.0 * c1) / gama;
				double term1 = pow((gama - 1) * vn1 + 2.0 * c1, 2);
				double term2 = 0.5 * (u1 * u1 + v1 * v1 - vn1 * vn1);
				//F+
				fluxVector11[IR] = fmass;
				fluxVector11[IU] = fmass * (u1 + nx * term0);
				fluxVector11[IV] = fmass * (v1 + ny * term0);
				//fluxVector11[IP] = fmass * h1;
				fluxVector11[IP] = fmass * (term1 / 2.0 / (pow(gama, 2) - 1) + term2);
			}

			if (m2 > 1.0)
			{
				fluxVector22[IR] = 0;
				fluxVector22[IU] = 0;
				fluxVector22[IV] = 0;
				fluxVector22[IP] = 0;
			}
			else if (m2 < -1.0)
			{
				Inviscid_Flux_F(fluxVector22, rho2, u2, v2, p2);
			}
			else
			{
				//F-
				double fmass = -rho2 * c2 * 0.25 * (m2 - 1) * (m2 - 1);
				double term0 = (-vn2 - 2.0 * c2) / gama;
				double term1 = pow((gama - 1) * vn2 - 2.0 * c2, 2);
				double term2 = 0.5 * (u2 * u2 + v2 * v2 - vn2 * vn2);

				fluxVector22[IR] = fmass;
				fluxVector22[IU] = fmass * (u2 + nx * term0);
				fluxVector22[IV] = fmass * (v2 + ny * term0);
				//fluxVector22[IP] = fmass * h2;
				fluxVector22[IP] = fmass * (term1 / 2.0 / (pow(gama, 2) - 1) + term2);
			}

			fluxVector[i][j][IR] = fluxVector11[IR] + fluxVector22[IR];
			fluxVector[i][j][IU] = fluxVector11[IU] + fluxVector22[IU];
			fluxVector[i][j][IV] = fluxVector11[IV] + fluxVector22[IV];
			fluxVector[i][j][IP] = fluxVector11[IP] + fluxVector22[IP];
		}
	}
#ifdef _OPENMP
}
#endif // _OPENMP
}

void Flux_Solver::VanLeer_Scheme_Y()
{
	VInt2D& marker = mesh->Get_Marker();

	double nx = 0.0, ny = 1.0;
#ifndef _OPENMP
	VDouble fluxVector11(num_of_prim_vars);
	VDouble fluxVector22(num_of_prim_vars);
	double rho1, u1, v1, p1;
	double rho2, u2, v2, p2;
#endif
#ifdef _OPENMP
#pragma omp parallel 
	{
		VDouble fluxVector11(num_of_prim_vars);
		VDouble fluxVector22(num_of_prim_vars);
		double rho1, u1, v1, p1;
		double rho2, u2, v2, p2;
#pragma omp	for
#endif
	for (int i = ist; i <= ied; i++)
	{
		for (int j = 0; j <= jed + 1; j++)//每个通量点的值都有了
		{
			if (marker[i][j] == 0) continue;

			ExtractValue(qField1[i][j], rho1, u1, v1, p1);

			double c1   = sqrt(fabs(gama * p1 / rho1));//声速取绝对值
			double vn1  = nx * u1 + ny * v1;
			double m1   = vn1 / c1;
			double h1   = Enthalpy(rho1, u1, v1, p1, gama);

			ExtractValue(qField2[i][j], rho2, u2, v2, p2);

			double c2   = sqrt(fabs(gama * p2 / rho2));//声速取绝对值
			double vn2  = nx * u2 + ny * v2;
			double m2   = vn2 / c2;
			double h2   = Enthalpy(rho2, u2, v2, p2, gama);

			if (m1 > 1.0)
			{
				Inviscid_Flux_G(fluxVector11, rho1, u1, v1, p1);
			}
			else if (m1 < -1.0)
			{
				fluxVector11[IR] = 0;
				fluxVector11[IU] = 0;
				fluxVector11[IV] = 0;
				fluxVector11[IP] = 0;
			}
			else
			{
				double fmass = rho1 * c1 * 0.25 * (m1 + 1) * (m1 + 1);
				double term0 = (-vn1 + 2.0 * c1) / gama;
				double term1 = pow((gama - 1) * vn1 + 2.0 * c1, 2);
				double term2 = 0.5 * (u1 * u1 + v1 * v1 - vn1 * vn1);
				//G+
				fluxVector11[IR] = fmass;
				fluxVector11[IU] = fmass * (u1 + nx * term0);
				fluxVector11[IV] = fmass * (v1 + ny * term0);
				//fluxVector11[IP] = fmass * h1;
				fluxVector11[IP] = fmass * (term1 / 2.0 / (pow(gama, 2) - 1) + term2);
			}

			if (m2 > 1.0)
			{
				fluxVector22[IR] = 0;
				fluxVector22[IU] = 0;
				fluxVector22[IV] = 0;
				fluxVector22[IP] = 0;
			}
			else if (m2 < -1.0)
			{
				Inviscid_Flux_G(fluxVector22, rho2, u2, v2, p2);
			}
			else
			{
				//G-
				double fmass = -rho2 * c2 * 0.25 * (m2 - 1) * (m2 - 1);
				double term0 = (-vn2 - 2.0 * c2) / gama;
				double term1 = pow((gama - 1) * vn2 - 2.0 * c2, 2);
				double term2 = 0.5 * (u2 * u2 + v2 * v2 - vn2 * vn2);

				fluxVector22[IR] = fmass;
				fluxVector22[IU] = fmass * (u2 + nx * term0);
				fluxVector22[IV] = fmass * (v2 + ny * term0);
				//fluxVector22[IP] = fmass * h2;
				fluxVector22[IP] = fmass * (term1 / 2.0 / (pow(gama, 2) - 1) + term2);
			}
			
			fluxVector[i][j][IR] = fluxVector11[IR] + fluxVector22[IR];
			fluxVector[i][j][IU] = fluxVector11[IU] + fluxVector22[IU];
			fluxVector[i][j][IV] = fluxVector11[IV] + fluxVector22[IV];
			fluxVector[i][j][IP] = fluxVector11[IP] + fluxVector22[IP];
		}
	}
#ifdef _OPENMP
	}
#endif // _OPENMP
}

void Flux_Solver::Steger_Warming_Scheme()
{
//	if (method_of_half_q == 1)		//插值之后计算半节点通量
//	{
		if (solve_direction == 'x')
		{
			this->Steger_Warming_Scheme_Interp_X();
		}
		else if (solve_direction == 'y')
		{
			this->Steger_Warming_Scheme_Interp_Y();
		}
//	}
}

void Flux_Solver::Steger_Warming_Scheme_X()
{
	double eps = 1e-4;
	VInt2D& marker = mesh->Get_Marker();

	//VDouble3D fluxVector1;
	//VDouble3D fluxVector2;
	//Allocate_3D_Vector(fluxVector1, total_points_x, total_points_y, num_of_prim_vars);
	//Allocate_3D_Vector(fluxVector2, total_points_x, total_points_y, num_of_prim_vars);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 0; i < ied + 2; i++)	//所有整数节点通量
		{
			if (marker[i][j] == 0) continue;

			double rho, u, v, p;
			ExtractValue(qField[i][j], rho, u, v, p);
			double a = sqrt(fabs(gama * p / rho));

			double lmd1 = u;
			double lmd3 = u - a;
			double lmd4 = u + a;

			VDouble lmd_m(3);//lamda-
			lmd_m[0] = 0.5 * (lmd1 - sqrt(lmd1 * lmd1 + eps * eps));
			lmd_m[1] = 0.5 * (lmd3 - sqrt(lmd3 * lmd3 + eps * eps));
			lmd_m[2] = 0.5 * (lmd4 - sqrt(lmd4 * lmd4 + eps * eps));

			Steger_Flux_F(fluxVector1[i][j], rho, u, v, p, lmd_m);

			VDouble lmd_p(3);//lamda+
			lmd_p[0] = 0.5 * (lmd1 + sqrt(lmd1 * lmd1 + eps * eps));
			lmd_p[1] = 0.5 * (lmd3 + sqrt(lmd3 * lmd3 + eps * eps));
			lmd_p[2] = 0.5 * (lmd4 + sqrt(lmd4 * lmd4 + eps * eps));

			Steger_Flux_F(fluxVector2[i][j], rho, u, v, p, lmd_p);
		}
	}
}

void Flux_Solver::Steger_Warming_Scheme_Y()
{
	double eps = 1e-4;
	VInt2D& marker = mesh->Get_Marker();
	//VDouble3D fluxVector1;
	//VDouble3D fluxVector2;
	//Allocate_3D_Vector(fluxVector1, total_points_x, total_points_y, num_of_prim_vars);
	//Allocate_3D_Vector(fluxVector2, total_points_x, total_points_y, num_of_prim_vars);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j = 0; j <= jed + 2; j++)
	{
		for (int i = ist; i <= ied; i++)
		{
			if (marker[i][j] == 0) continue;

			double rho, u, v, p;
			ExtractValue(qField[i][j], rho, u, v, p);

			double a = sqrt(fabs(gama * p / rho));

			double mu1 = v;
			double mu3 = v - a;
			double mu4 = v + a;

			VDouble mu_m(3);//mu-
			mu_m[0] = 0.5 * (mu1 - sqrt(mu1 * mu1 + eps * eps));
			mu_m[1] = 0.5 * (mu3 - sqrt(mu3 * mu3 + eps * eps));
			mu_m[2] = 0.5 * (mu4 - sqrt(mu4 * mu4 + eps * eps));

			Steger_Flux_G(fluxVector1[i][j], rho, u, v, p, mu_m);

			VDouble mu_p(3);//mu+
			mu_p[0] = 0.5 * (mu1 + sqrt(mu1 * mu1 + eps * eps));
			mu_p[1] = 0.5 * (mu3 + sqrt(mu3 * mu3 + eps * eps));
			mu_p[2] = 0.5 * (mu4 + sqrt(mu4 * mu4 + eps * eps));

			Steger_Flux_G(fluxVector2[i][j], rho, u, v, p, mu_p);
		}
	}
}

void Flux_Solver::Steger_Warming_Scheme_Interp_X()
{
	double eps = 1e-4;
	VInt2D& marker = mesh->Get_Marker();
#ifndef _OPENMP
	VDouble fluxVector11(num_of_prim_vars);
	VDouble fluxVector22(num_of_prim_vars);
#endif
#ifdef _OPENMP
#pragma omp parallel 
{
		VDouble fluxVector11(num_of_prim_vars);
		VDouble fluxVector22(num_of_prim_vars);
#pragma omp	for
#endif
	for (int j = jst; j <= jed; j++)
	{
		for (int i = 0; i <= ied + 1; i++)//每个通量点的值都有了
		{
			if (marker[i][j] == 0) continue;

			double rho, u, v, p;
			ExtractValue(qField1[i][j], rho, u, v, p);

			double a    = sqrt(fabs(gama * p / rho));//声速取绝对值

			double lmd1 = u;
			double lmd3 = u - a;
			double lmd4 = u + a;

			VDouble lmd_p(3);//q1对应于lamda+
			lmd_p[0] = 0.5 * (lmd1 + sqrt(lmd1 * lmd1 + eps * eps));
			lmd_p[1] = 0.5 * (lmd3 + sqrt(lmd3 * lmd3 + eps * eps));
			lmd_p[2] = 0.5 * (lmd4 + sqrt(lmd4 * lmd4 + eps * eps));
			
			Steger_Flux_F(fluxVector11, rho, u, v, p, lmd_p);	//F+

			ExtractValue(qField2[i][j], rho, u, v, p);

			a    = sqrt(fabs(gama * p / rho));//声速取绝对值

			lmd1 = u;
			lmd3 = u - a;
			lmd4 = u + a;

			VDouble lmd_m(3);//q2对应于lamda-
			lmd_m[0] = 0.5 * (lmd1 - sqrt(lmd1 * lmd1 + eps * eps));
			lmd_m[1] = 0.5 * (lmd3 - sqrt(lmd3 * lmd3 + eps * eps));
			lmd_m[2] = 0.5 * (lmd4 - sqrt(lmd4 * lmd4 + eps * eps));
			
			Steger_Flux_F(fluxVector22, rho, u, v, p, lmd_m);	//F-

			fluxVector[i][j][IR] = fluxVector11[IR] + fluxVector22[IR];
			fluxVector[i][j][IU] = fluxVector11[IU] + fluxVector22[IU];
			fluxVector[i][j][IV] = fluxVector11[IV] + fluxVector22[IV];
			fluxVector[i][j][IP] = fluxVector11[IP] + fluxVector22[IP];
		}
	}
#ifdef _OPENMP
}
#endif // _OPENMP
}

void Flux_Solver::Steger_Warming_Scheme_Interp_Y()
{
	double eps = 1e-4;
	VInt2D& marker = mesh->Get_Marker();
#ifndef _OPENMP
	VDouble fluxVector11(num_of_prim_vars);
	VDouble fluxVector22(num_of_prim_vars);
#endif
#ifdef _OPENMP
#pragma omp parallel 
	{
		VDouble fluxVector11(num_of_prim_vars);
		VDouble fluxVector22(num_of_prim_vars);
#pragma omp	for
#endif
	for (int j = 0; j <= jed + 1; j++)  //每个通量点的值都有了
	{
		for (int i = ist; i <= ied; i++)
		{
			if (marker[i][j] == 0) continue;

			double rho, u, v, p;
			ExtractValue(qField1[i][j], rho, u, v, p);

			double a   = sqrt(fabs(gama * p / rho));//声速取绝对值

			double mu1 = v;
			double mu3 = v - a;
			double mu4 = v + a;

			VDouble mu_p(3);//q1对应于mu+
			mu_p[0] = 0.5 * (mu1 + sqrt(mu1 * mu1 + eps * eps));
			mu_p[1] = 0.5 * (mu3 + sqrt(mu3 * mu3 + eps * eps));
			mu_p[2] = 0.5 * (mu4 + sqrt(mu4 * mu4 + eps * eps));
			
			Steger_Flux_G(fluxVector11, rho, u, v, p, mu_p);	//G+

			ExtractValue(qField2[i][j], rho, u, v, p);

			a   = sqrt(fabs(gama * p / rho));//声速取绝对值

			mu1 = v;
			mu3 = v - a;
			mu4 = v + a;

			VDouble mu_m(3);//q2对应于mu-
			mu_m[0] = 0.5 * (mu1 - sqrt(mu1 * mu1 + eps * eps));
			mu_m[1] = 0.5 * (mu3 - sqrt(mu3 * mu3 + eps * eps));
			mu_m[2] = 0.5 * (mu4 - sqrt(mu4 * mu4 + eps * eps));
			
			Steger_Flux_G(fluxVector22, rho, u, v, p, mu_m);	//G-

			fluxVector[i][j][IR] = fluxVector11[IR] + fluxVector22[IR];
			fluxVector[i][j][IU] = fluxVector11[IU] + fluxVector22[IU];
			fluxVector[i][j][IV] = fluxVector11[IV] + fluxVector22[IV];
			fluxVector[i][j][IP] = fluxVector11[IP] + fluxVector22[IP];
		}
	}
#ifdef _OPENMP
	}
#endif // _OPENMP
}

void Flux_Solver::Steger_Flux_F(VDouble& fluxVector, double rho, double u, double v, double p, VDouble lmd)
{
	double a = sqrt(fabs(gama * p / rho));//声速取绝对值
	double h = Enthalpy(rho, u, v, p, gama);

	fluxVector[IR] = rho / (2 * gama) * (2 * (gama - 1) *	  lmd[0] +			 lmd[1] +			    lmd[2]);
	fluxVector[IU] = rho / (2 * gama) * (2 * (gama - 1) * u * lmd[0] + (u - a) * lmd[1] + (u + a)	  * lmd[2]);
	fluxVector[IV] = rho / (2 * gama) * (2 * (gama - 1) * v * lmd[0] +		 v * lmd[1] +		    v * lmd[2]);
	fluxVector[IP] = rho / (2 * gama) * ((gama - 1) *		   (u * u + v * v) * lmd[0] + (h - a * u) * lmd[1] + (h + a * u) * lmd[2]);
}	

void Flux_Solver::Steger_Flux_G(VDouble& fluxVector, double rho, double u, double v, double p, VDouble mu)
{
	double a = sqrt(fabs(gama * p / rho));//声速取绝对值
	double h = Enthalpy(rho, u, v, p, gama);

	fluxVector[IR] = rho / (2 * gama) * (2 * (gama - 1) *	  mu[0] +			mu[1] +				  mu[2]);
	fluxVector[IU] = rho / (2 * gama) * (2 * (gama - 1) * u * mu[0] +		u * mu[1] +			  u * mu[2]);
	fluxVector[IV] = rho / (2 * gama) * (2 * (gama - 1) * v * mu[0] + (v - a) * mu[1] +		(v + a) * mu[2]);
	fluxVector[IP] = rho / (2 * gama) * (    (gama - 1) *     (v * v + u * u) * mu[0] + (h - a * v) * mu[1] + (h + a * v) * mu[2]);
}

void Flux_Solver::Inviscid_Flux_F(VDouble& fluxVector, double rho, double u, double v, double p)
{
	double E = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
	fluxVector[IR] = rho * u;
	fluxVector[IU] = rho * u * u + p;
	fluxVector[IV] = rho * u * v;
	fluxVector[IP] = (E + p) * u;
}

void Flux_Solver::Inviscid_Flux_G(VDouble& fluxVector, double rho, double u, double v, double p)
{
	double E = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
	fluxVector[IR] = rho * v;
	fluxVector[IU] = rho * u * v;
	fluxVector[IV] = rho * v * v + p;
	fluxVector[IP] = (E + p) * v;
}

double Flux_Solver::Enthalpy(double rho, double u, double v, double p, double gama)
{
	double E = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
	return (E + p) / rho;
}

void Flux_Solver::EntropyFix(double& lamda1, double& lamda2, double& lamda3, double& lamda4)
{
	this->EntropyFix_Harten(lamda1);
	this->EntropyFix_Harten(lamda2);
	this->EntropyFix_Harten(lamda3);
	this->EntropyFix_Harten(lamda4);
}

void Flux_Solver::EntropyFix_Harten(double &lamda)
{
	double eps = entropy_fix_coeff;
	lamda = fabs(lamda) > eps ? fabs(lamda) : ((lamda * lamda + eps * eps) / 2.0 / eps);
}

void Flux_Solver::Compute_Jacobian(VDouble2D& Jacobian, double u, double v, double a, double h)
{
	if (solve_direction == 'x')
	{
		double lamda1 = u;
		double lamda2 = u;
		double lamda3 = u - a;
		double lamda4 = u + a;

		EntropyFix(lamda1, lamda2, lamda3, lamda4);

		Compute_Jacobian_X(Jacobian, u, v, a, h, lamda1, lamda2, lamda3, lamda4);
	}
	else if (solve_direction == 'y')
	{
		double lamda1 = v;
		double lamda2 = v;
		double lamda3 = v - a;
		double lamda4 = v + a;

		EntropyFix(lamda1, lamda2, lamda3, lamda4);

		Compute_Jacobian_Y(Jacobian, u, v, a, h, lamda1, lamda2, lamda3, lamda4);
	}	
}

void Flux_Solver::Compute_Jacobian_X(VDouble2D& Jacobian, double u, double v, double a, double h,
	double lamda1, double lamda2, double lamda3, double lamda4)
{
	VDouble2D Lmdx = { {lamda1,0,0,0},{0,lamda2,0,0},{0,0,lamda3,0},{0,0,0,lamda4} };

	VDouble2D Rx = { {1,						 0,   1,			1			},
					 {u,						 0,   u - a,		u + a		},
					 {0,						 1,   v,			v,			},
					 {(u * u - v * v) / 2,     v,     h - a * u,    h + a * u	} };
	
	double b2 = (gama - 1) / a / a;
	double b1 = b2 * (u * u + v * v) / 2.0;

	VDouble2D Lx = { {1 - b1,			   b2 * u,				  b2 * v,		   - b2		},
					 {-b1 * v,			   b2 * u * v,            1 + b2 * v * v,  - b2 * v	},
					 {(b1 + u / a) / 2, - (b2 * u + 1 / a) / 2, - b2 * v / 2,        b2 / 2	},
					 {(b1 - u / a) / 2, - (b2 * u - 1 / a) / 2, - b2 * v / 2,		 b2 / 2	}};

	VDouble2D tmp1;
	Allocate_2D_Vector(tmp1, 4, 4);

	MatrixMultiply(Rx,   Lmdx, tmp1,     4, 4, 4);
	MatrixMultiply(tmp1, Lx,   Jacobian, 4, 4, 4);
}

void Flux_Solver::Compute_Jacobian_Y(VDouble2D& Jacobian, double u, double v, double a, double h,
	double lamda1, double lamda2, double lamda3, double lamda4)
{
	VDouble2D Lmdy = { {lamda1,0,0,0},{0,lamda2,0,0},{0,0,lamda3,0},{0,0,0,lamda4} };

	VDouble2D Ry = { {0,						 1,     1,			1		  },
					 {1,						 0,     u,			u  		  },
					 {0,						 v,     v - a,		v + a,	  },
					 {u,       (v * v - u * u) / 2,     h - a * v,  h + a * v } };

	double b2 = (gama - 1) / a / a;
	double b1 = b2 * (u * u + v * v) / 2.0;

	VDouble2D Ly = { {         -b1 * u,		1 + b2 * u * u,               b2 * u * v,    -b2 * u	},
					 {          1 - b1,		        b2 * u,				      b2 * v,	 -b2		},
					 {(b1 + v / a) / 2,		   -b2 * u / 2,    -(b2 * v + 1 / a) / 2,     b2 / 2	},
					 {(b1 - v / a) / 2,        -b2 * u / 2,    -(b2 * v - 1 / a) / 2,	  b2 / 2	} };

	VDouble2D tmp1;
	Allocate_2D_Vector(tmp1, 4, 4);

	MatrixMultiply(Ry, Lmdy, tmp1, 4, 4, 4);
	MatrixMultiply(tmp1, Ly, Jacobian, 4, 4, 4);
}

