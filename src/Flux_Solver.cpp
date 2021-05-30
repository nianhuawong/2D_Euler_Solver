#include <iostream>
#include "Global.h"
#include "QlQr_Solver.h"
#include "Flux_Solver.h"
#include "2D_Euler_Solver.h"
#include "Geometry.h"

vector< vector< vector< double > > > fluxVector;

void Solve_Flux()
{
	auto* fluxVector = new Flux_Solver();

	fluxVector->Solve_Flux();

	delete fluxVector;
}

Flux_Solver::Flux_Solver()
{
	this->gama = 1.4;

	fluxVector.resize(num_half_point_x);
	this->fluxVector1.resize (num_half_point_x);
	this->fluxVector2.resize (num_half_point_x);
	for (int i = 0; i < num_half_point_x; i++)
	{
		Allocate_2D_Vector(fluxVector [i], num_half_point_y, num_of_prim_vars);
		Allocate_2D_Vector(fluxVector1[i], num_half_point_y, num_of_prim_vars);
		Allocate_2D_Vector(fluxVector2[i], num_half_point_y, num_of_prim_vars);
	}

	Allocate_2D_Vector(Jacobian_A, num_of_prim_vars, num_of_prim_vars);
}

void Flux_Solver::Solve_Flux()
{
	if (method_of_flux == 1)
	{
		Flux_LR_Roe();
		Roe_Scheme();
	}
	else if (method_of_flux == 2)
	{
		//WENO+Steger-Warming
	}
	else if (method_of_flux == 3)
	{
		//WCNS+Roe
	}
	else
	{
		cout << "flux计算方法选择错误，请检查！" << endl;
	}

}
void Flux_Solver::Flux_LR_Roe()
{
	if (solve_direction == 'x')
	{
		this->Flux_LR_Roe_X();
	}
	else if (solve_direction == 'y')
	{
		this->Flux_LR_Roe_Y();
	}
	else
	{
		cout << "出错，请检查！" << endl;
	}
}

void Flux_Solver::Flux_LR_Roe_X()
{
	for (int j = 0; j < num_half_point_y; j++)
	{
		for (int i = 0; i < num_half_point_x; i++)
		{
			double rho1 = qField1[IR][i][j];
			double u1   = qField1[IU][i][j];
			double v1	= qField1[IV][i][j];
			double E1   = qField1[IP][i][j];
			double p1	= Energy_2_Pressure(E1, rho1, u1, v1);

			Inviscid_Flux_F(fluxVector1[i][j], rho1, u1, v1, p1);

			double rho2 = qField2[IR][i][j];
			double u2   = qField2[IU][i][j];
			double v2   = qField2[IV][i][j];
			double E2   = qField2[IP][i][j];
			double p2	= Energy_2_Pressure(E2, rho2, u2, v2);

			Inviscid_Flux_F(fluxVector2[i][j], rho2, u2, v2, p2);
		}
	}
}

void Flux_Solver::Flux_LR_Roe_Y()
{
	for (int i = 0; i < num_half_point_x; i++)
	{
		for (int j = 0; j < num_half_point_y; j++)
		{
			double rho1 = qField1[IR][i][j];
			double u1 = qField1[IU][i][j];
			double v1 = qField1[IV][i][j];
			double E1 = qField1[IP][i][j];
			double p1 = Energy_2_Pressure(E1, rho1, u1, v1);

			Inviscid_Flux_G(fluxVector1[i][j], rho1, u1, v1, p1);

			double rho2 = qField2[IR][i][j];
			double u2 = qField2[IU][i][j];
			double v2 = qField2[IV][i][j];
			double E2 = qField2[IP][i][j];
			double p2 = Energy_2_Pressure(E2, rho2, u2, v2);

			Inviscid_Flux_G(fluxVector2[i][j], rho2, u2, v2, p2);
		}
	}
}

void Flux_Solver::Flux_LR_Steger_Warming()
{

}

void Flux_Solver::Inviscid_Flux_F(vector< double >& fluxVector, double rho, double u, double v, double p)
{
	double E = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
	fluxVector[IR] = rho * u;
	fluxVector[IU] = rho * u * u + p;
	fluxVector[IV] = rho * u * v;
	fluxVector[IP] = (E + p) * u;
}

void Flux_Solver::Inviscid_Flux_G(vector< double >& fluxVector, double rho, double u, double v, double p)
{
	double E = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
	fluxVector[IR] = rho * v;
	fluxVector[IU] = rho * u * v;
	fluxVector[IV] = rho * v * v + p;
	fluxVector[IP] = (E + p) * v;
}

void Flux_Solver::Roe_Scheme()
{
	for (int j = 0; j < num_half_point_y; j++)
	{
		for (int i = 0; i < num_half_point_x; i++)
		{
			double rho1 = qField1[IR][i][j];
			double u1   = qField1[IU][i][j];
			double v1	= qField1[IV][i][j];
			double E1   = qField1[IP][i][j];
			double p1   = Energy_2_Pressure(E1, rho1, u1, v1);
			double H1   = Enthalpy(rho1, p1, gama);

			double rho2 = qField2[IR][i][j];
			double u2   = qField2[IU][i][j];
			double v2   = qField2[IV][i][j];
			double E2   = qField2[IP][i][j];
			double p2   = Energy_2_Pressure(E2, rho2, u2, v2);
			double H2   = Enthalpy(rho2, p2, gama);

			double D = sqrt(rho2 / rho1);

			double rho_roe = sqrt(rho1 * rho2);
			double u_roe = (u1 + u2 * D) / (1 + D);
			double v_roe = (v1 + v2 * D) / (1 + D);
			double H_roe = (H1 + H2 * D) / (1 + D);
			double c_roe = sqrt((gama - 1) * (H_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe)));

			double lamda1 = u_roe;
			double lamda2 = u_roe;
			double lamda3 = u_roe - c_roe;
			double lamda4 = u_roe + c_roe;

			EntropyFix(&lamda1, &lamda2, &lamda3, &lamda4);	

			Compute_Jacobian(Jacobian_A, u_roe, v_roe, c_roe, H_roe, lamda1, lamda2, lamda3, lamda4);

			vector< double > dq(num_of_prim_vars);
			dq[IR] = qField2[IR][i][j] - qField1[IR][i][j];
			dq[IU] = qField2[IU][i][j] - qField1[IU][i][j];
			dq[IV] = qField2[IV][i][j] - qField1[IV][i][j];
			dq[IP] = qField2[IP][i][j] - qField1[IP][i][j];

			vector<double> flux_add(num_of_prim_vars);
			MatrixMultiply(Jacobian_A, dq, flux_add, 4, 4);

			fluxVector[i][j][IR] = 0.5 * (fluxVector1[i][j][IR] + fluxVector2[i][j][IR] + flux_add[IR]);
			fluxVector[i][j][IU] = 0.5 * (fluxVector1[i][j][IU] + fluxVector2[i][j][IU] + flux_add[IU]);
			fluxVector[i][j][IV] = 0.5 * (fluxVector1[i][j][IV] + fluxVector2[i][j][IV] + flux_add[IV]);
			fluxVector[i][j][IP] = 0.5 * (fluxVector1[i][j][IP] + fluxVector2[i][j][IP] + flux_add[IP]);
		}
	}
}

double Flux_Solver::Enthalpy(double rho, double p, double gama)
{
	double e = p / (gama - 1) / rho;
	return e + p / rho;
}

void Flux_Solver::EntropyFix(double* lamda1, double* lamda2, double* lamda3, double* lamda4)
{
	*lamda1 = this->EntropyFix_Harten(*lamda1);
	*lamda2 = this->EntropyFix_Harten(*lamda2);
	*lamda3 = this->EntropyFix_Harten(*lamda3);
	*lamda4 = this->EntropyFix_Harten(*lamda4);
}

double Flux_Solver::EntropyFix_Harten(double lamda)
{
	double eps = entropy_fix_coeff;
	return abs(lamda) > eps ? abs(lamda) : ((lamda * lamda + eps * eps) / 2.0 / eps);
}

void Flux_Solver::Compute_Jacobian(vector < vector<double> >& Jacobian, double u, double v, double a, double h,
	double lamda1, double lamda2, double lamda3, double lamda4)
{
	vector < vector< double > > Lmdx = { {lamda1,0,0,0},{0,lamda2,0,0},{0,0,lamda3,0},{0,0,0,lamda4} };

	vector < vector< double > > Rx = { {1,						 0,   1,			1			},
									   {u,						 0,   u - a,		u + a		},
									   {0,						 1,   v,			v,			},
									   {(u * u - v * v) / 2,     v,   h - a * u,    h + a * u	} };
	double b2 = (gama - 1) / a / a;
	double b1 = b2 * (u * u + v * v) / 2.0;

	vector < vector< double > > Lx = { {1 - b1,				 b2 * u,				b2 * v,			 - b2		},
									   {-b1 * v,			 b2 * u * v,            1 + b2 * v * v,  - b2 * v	},
									   {(b1 + u / a) / 2, - (b2 * u + 1 / a) / 2, - b2 * v / 2,      b2 / 2		},
									   {(b1 - u / a) / 2, - (b2 * u - 1 / a) / 2, - b2 * v / 2,		 b2 / 2		}};

	vector < vector< double > > tmp1;
	Allocate_2D_Vector(tmp1, 4, 4);

	MatrixMultiply(Rx,   Lmdx, tmp1,     4, 4, 4);
	MatrixMultiply(tmp1, Lx,   Jacobian, 4, 4, 4);
}

