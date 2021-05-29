#pragma once
#include<vector>
using namespace std;

class Half_Node_Flux_Solver
{
public:
	Half_Node_Flux_Solver();
	~Half_Node_Flux_Solver() {};
protected:
	double gama;

	vector< vector< vector< double > > > half_node_flux_l;
	vector< vector< vector< double > > > half_node_flux_r;
	vector < vector<double> > Jacobian_A;

public:
	void Half_Node_Flux();
	void Half_Node_Flux_LR_Roe();
	void Half_Node_Flux_LR_Steger_Warming();

protected:
	void Inviscid_Flux_F(vector< double >& fluxVector, double rho, double u, double v, double p);
	void Inviscid_Flux_G(vector< double >& fluxVector, double rho, double u, double v, double p);
	void Roe_Scheme();

	double Enthalpy(double rho, double p, double gama);
	void EntropyFix(double &lamda1, double& lamda2, double& lamda3, double& lamda4);
	double EntropyFix_Harten(double lamda);

	void Compute_Jacobian(vector < vector<double> > &Jacobian, double u,  double v,  double c,  double H, 
															   double lamda1, double lamda2, double lamda3, double lamda4);
};
void Half_Node_Flux();