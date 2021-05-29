#include <iostream>
#include "Global.h"
#include "QlQr_Solver.h"
#include "Flux_Solver.h"
#include "2D_Euler_Solver.h"
#include "Geometry.h"

using namespace GLOBAL;

void GLOBAL::Solve_Flux()
{
	auto* fluxVector = new Flux_Solver();

	fluxVector->Solve_Flux();

	delete fluxVector;
}

Flux_Solver::Flux_Solver()
{
	this->gama = 1.4;

	GLOBAL::fluxVector.resize(num_half_point_x);
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
	Flux_LR_Roe();

	Roe_Scheme();
}

void Flux_Solver::Flux_LR_Roe()
{
	for (int j = 0; j < num_half_point_y; j++)
	{
		for (int i = 0; i < num_half_point_x; i++)
		{
			double rho1 = qField1[IR][i][j];
			double u1   = qField1[IU][i][j];
			double v1	= qField1[IV][i][j];
			double p1   = qField1[IP][i][j];

			Inviscid_Flux_F(fluxVector1[i][j], rho1, u1, v1, p1);

			double rho2 = qField2[IR][i][j];
			double u2   = qField2[IU][i][j];
			double v2   = qField2[IV][i][j];
			double p2   = qField2[IP][i][j];
			Inviscid_Flux_F(fluxVector2[i][j], rho2, u2, v2, p2);
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
			double p1   = qField1[IP][i][j];		
			double H1   = Enthalpy(rho1, p1, gama);

			double rho2 = qField2[IR][i][j];
			double u2   = qField2[IU][i][j];
			double v2   = qField2[IV][i][j];
			double p2   = qField2[IP][i][j];
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

			vector< double >dq(num_of_prim_vars);
			dq[IR] = qField2[IR][i][j] - qField1[IR][i][j];
			dq[IU] = qField2[IU][i][j] - qField1[IU][i][j];
			dq[IV] = qField2[IV][i][j] - qField1[IV][i][j];
			dq[IP] = qField2[IP][i][j] - qField1[IP][i][j];

			vector<double> flux_add(num_of_prim_vars);
			flux_add[IR] = Jacobian_A[IR][0] * dq[0] + Jacobian_A[IR][1] * dq[1] + Jacobian_A[IR][2] * dq[2] + Jacobian_A[IR][3] * dq[3];
			flux_add[IU] = Jacobian_A[IU][0] * dq[0] + Jacobian_A[IU][1] * dq[1] + Jacobian_A[IU][2] * dq[2] + Jacobian_A[IU][3] * dq[3];
			flux_add[IV] = Jacobian_A[IV][0] * dq[0] + Jacobian_A[IV][1] * dq[1] + Jacobian_A[IV][2] * dq[2] + Jacobian_A[IV][3] * dq[3];
			flux_add[IP] = Jacobian_A[IP][0] * dq[0] + Jacobian_A[IP][1] * dq[1] + Jacobian_A[IP][2] * dq[2] + Jacobian_A[IP][3] * dq[3];

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
	//double gm = gama;
	//Jacobian[0][0] = lamda3 * (u / (2 * a) + ((u * u + v * v) * (gm - 1)) / (4 * a * a)) - lamda1 * (((u * u + v * v) * (gm - 1)) / (2 * a * a) - 1) - lamda4 * (u / (2 * a) - ((u * u + v * v) * (gm - 1)) / (4 * a * a));
	//Jacobian[0][1] = lamda4 * (1 / (2 * a) - (u * (gm - 1)) / (2 * a * a)) - lamda3 * (1 / (2 * a) + (u * (gm - 1)) / (2 * a * a)) + (lamda1 * u * (gm - 1)) / a / a;
	//Jacobian[0][2] = 
}