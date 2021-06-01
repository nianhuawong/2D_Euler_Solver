#include <iostream>
#include "Global.h"
#include "QlQr_Solver.h"
#include "Flux_Solver.h"
#include "2D_Euler_Solver.h"
#include "Geometry.h"

VDouble3D fluxVector;

void Solve_Flux()
{
	auto* fluxVector = new Flux_Solver();

	fluxVector->Solve_Flux();

	delete fluxVector;
}

Flux_Solver::Flux_Solver()
{
	this->gama = 1.4;

	Get_IJK_Region(ist, ied, jst, jed);

	Allocate_3D_Vector(fluxVector,  num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(fluxVector1, num_half_point_x, num_half_point_y, num_of_prim_vars);
	Allocate_3D_Vector(fluxVector2, num_half_point_x, num_half_point_y, num_of_prim_vars);

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
		Flux_LR_Steger_Warming();
		Steger_Warming_Scheme();	//二阶Steger_Warming
	}
	else if (method_of_flux == 3)
	{
		//WENO+Steger-Warming
		Flux_LR_Steger_Warming(); 
		//WENO通量插值如下：

	}
	else if (method_of_flux == 4)
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
	VInt2D& marker = mesh->Get_Marker();
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;

			double rho1 = qField1[i][j][IR];
			double u1   = qField1[i][j][IU];
			double v1	= qField1[i][j][IV];
			double p1   = qField1[i][j][IP];

			Inviscid_Flux_F(fluxVector1[i][j], rho1, u1, v1, p1);

			double rho2 = qField2[i][j][IR];
			double u2   = qField2[i][j][IU];
			double v2   = qField2[i][j][IV];
			double p2   = qField2[i][j][IP];

			Inviscid_Flux_F(fluxVector2[i][j], rho2, u2, v2, p2);
		}
	}
}

void Flux_Solver::Flux_LR_Roe_Y()
{
	VInt2D& marker = mesh->Get_Marker();
	for (int i = ist; i < ied - 1; i++)
	{
		for (int j = jst; j < jed - 1; j++)
		{
			if (marker[i][j] == 0) continue;

			double rho1 = qField1[i][j][IR];
			double u1   = qField1[i][j][IU];
			double v1   = qField1[i][j][IV];
			double p1   = qField1[i][j][IP];

			Inviscid_Flux_G(fluxVector1[i][j], rho1, u1, v1, p1);

			double rho2 = qField2[i][j][IR];
			double u2   = qField2[i][j][IU];
			double v2   = qField2[i][j][IV];
			double p2   = qField2[i][j][IP];

			Inviscid_Flux_G(fluxVector2[i][j], rho2, u2, v2, p2);
		}
	}
}

void Flux_Solver::Flux_LR_Steger_Warming()
{
	if (solve_direction == 'x')
	{
		this->Flux_LR_Steger_Warming_X();
	}
	else if (solve_direction == 'y')
	{
		this->Flux_LR_Steger_Warming_Y();
	}
	else
	{
		cout << "出错，请检查！" << endl;
	}
}

void Flux_Solver::Steger_Warming_Scheme()
{
	VInt2D& marker = mesh->Get_Marker();
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;

			fluxVector[i][j][IR] = fluxVector1[i][j][IR] + fluxVector2[i][j][IR];
			fluxVector[i][j][IU] = fluxVector1[i][j][IU] + fluxVector2[i][j][IU];
			fluxVector[i][j][IV] = fluxVector1[i][j][IV] + fluxVector2[i][j][IV];
			fluxVector[i][j][IP] = fluxVector1[i][j][IP] + fluxVector2[i][j][IP];
		}
	}
}

void Flux_Solver::Flux_LR_Steger_Warming_X()
{
	double eps = 1e-4;
	VInt2D& marker = mesh->Get_Marker();
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;

			double rho = qField[i][j][IR];
			double u = qField[i][j][IU];
			double v = qField[i][j][IV];
			double p = qField[i][j][IP];
			double a = sqrt(gama * p / rho);

			double lmd1 = u;
			double lmd3 = u - a;
			double lmd4 = u + a;

			VDouble lmd_p(3);//lamda+
			lmd_p[0] = 0.5 * (lmd1 + sqrt(lmd1 * lmd1 + eps * eps));
			lmd_p[1] = 0.5 * (lmd3 + sqrt(lmd3 * lmd3 + eps * eps));
			lmd_p[2] = 0.5 * (lmd4 + sqrt(lmd4 * lmd4 + eps * eps));

			VDouble lmd_m(3);//lamda-
			lmd_p[0] = 0.5 * (lmd1 - sqrt(lmd1 * lmd1 + eps * eps));
			lmd_m[1] = 0.5 * (lmd3 - sqrt(lmd3 * lmd3 + eps * eps));
			lmd_m[2] = 0.5 * (lmd4 - sqrt(lmd4 * lmd4 + eps * eps));

			Steger_Flux_F(fluxVector1[i][j], rho, u, v, p, lmd_p);
			Steger_Flux_F(fluxVector2[i][j], rho, u, v, p, lmd_m);
		}
	}
}

void Flux_Solver::Flux_LR_Steger_Warming_Y()
{
	double eps = 1e-4;
	VInt2D& marker = mesh->Get_Marker();
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;

			double rho = qField[i][j][IR];
			double u = qField[i][j][IU];
			double v = qField[i][j][IV];
			double p = qField[i][j][IP];
			double a = sqrt(gama * p / rho);

			double mu1 = v;
			double mu3 = v - a;
			double mu4 = v + a;

			VDouble mu_p(3);//mu+
			mu_p[0] = 0.5 * (mu1 + sqrt(mu1 * mu1 + eps * eps));
			mu_p[1] = 0.5 * (mu3 + sqrt(mu3 * mu3 + eps * eps));
			mu_p[2] = 0.5 * (mu4 + sqrt(mu4 * mu4 + eps * eps));

			VDouble mu_m(3);//mu-
			mu_m[0] = 0.5 * (mu1 - sqrt(mu1 * mu1 + eps * eps));
			mu_m[1] = 0.5 * (mu3 - sqrt(mu3 * mu3 + eps * eps));
			mu_m[2] = 0.5 * (mu4 - sqrt(mu4 * mu4 + eps * eps));

			Steger_Flux_G(fluxVector1[i][j], rho, u, v, p, mu_p);
			Steger_Flux_G(fluxVector2[i][j], rho, u, v, p, mu_m);
		}
	}
}

void Flux_Solver::Steger_Flux_F(VDouble& fluxVector, double rho, double u, double v, double p, VDouble lmd)
{
	double a = sqrt(gama * p / rho);
	double h = Enthalpy(rho, u, v, p, gama);

	fluxVector[IR] = rho / (2 * gama) * (2 * (gama - 1) * lmd[0] + lmd[2] + lmd[3]);
	fluxVector[IU] = rho / (2 * gama) * (2 * (gama - 1) * u * lmd[0] + (u - a) * lmd[2] + (u + a) * lmd[3]);
	fluxVector[IV] = rho / (2 * gama) * (2 * (gama - 1) * v * lmd[0] + v * lmd[2] + v * lmd[3]);
	fluxVector[IP] = rho / (2 * gama) * ((gama - 1) * (pow(u, 2) + pow(v, 2)) * lmd[0] + (h - a * u) * lmd[2] + (h + a * u) * lmd[3]);
}	

void Flux_Solver::Steger_Flux_G(VDouble& fluxVector, double rho, double u, double v, double p, VDouble mu)
{
	double a = sqrt(gama * p / rho);
	double h = Enthalpy(rho, u, v, p, gama);

	fluxVector[IR] = rho / (2 * gama) * (2 * (gama - 1) * mu[0] + mu[2] + mu[3]);
	fluxVector[IU] = rho / (2 * gama) * (2 * (gama - 1) * u * mu[0] + u * mu[2] + u * mu[3]);
	fluxVector[IV] = rho / (2 * gama) * (2 * (gama - 1) * v * mu[0] + (v - a) * mu[2] + (v + a) * mu[3]);
	fluxVector[IP] = rho / (2 * gama) * ((gama - 1) * (pow(u, 2) + pow(v, 2)) * mu[0] + (h - a * v) * mu[2] + (h + a * v) * mu[3]);
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

void Flux_Solver::Roe_Scheme()
{
	VInt2D& marker = mesh->Get_Marker();
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;

			double rho1 = qField1[i][j][IR];
			double u1   = qField1[i][j][IU];
			double v1	= qField1[i][j][IV];
			double p1   = qField1[i][j][IP];
			double H1  = Enthalpy(rho1, u1, v1, p1, gama);

			double rho2 = qField2[i][j][IR];
			double u2   = qField2[i][j][IU];
			double v2   = qField2[i][j][IV];
			double p2   = qField2[i][j][IP];
			double H2  = Enthalpy(rho2, u2, v2, p2, gama);

			double D = sqrt(rho2 / rho1);

			double rho_roe = sqrt(rho1 * rho2);
			double u_roe = (u1 + u2 * D) / (1 + D);
			double v_roe = (v1 + v2 * D) / (1 + D);
			double H_roe = (H1 + H2 * D) / (1 + D);
			double c2_roe = (gama - 1) * (H_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe));
			double c_roe = sqrt(abs(c2_roe));

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

			fluxVector[i][j][IR] = 0.5 * (fluxVector1[i][j][IR] + fluxVector2[i][j][IR] - flux_add[IR]);
			fluxVector[i][j][IU] = 0.5 * (fluxVector1[i][j][IU] + fluxVector2[i][j][IU] - flux_add[IU]);
			fluxVector[i][j][IV] = 0.5 * (fluxVector1[i][j][IV] + fluxVector2[i][j][IV] - flux_add[IV]);
			fluxVector[i][j][IP] = 0.5 * (fluxVector1[i][j][IP] + fluxVector2[i][j][IP] - flux_add[IP]);
		}
	}
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
	lamda = abs(lamda) > eps ? abs(lamda) : ((lamda * lamda + eps * eps) / 2.0 / eps);
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

