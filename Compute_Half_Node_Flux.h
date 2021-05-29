#pragma once
#include<vector>
using namespace std;

class Half_Node_Flux_Solver
{
public:
	Half_Node_Flux_Solver();
	~Half_Node_Flux_Solver() {};
protected:
	int M_Dim, N_Dim;
	int method_of_flux;

	vector< vector< vector< double > > > half_node_flux_l;
	vector< vector< vector< double > > > half_node_flux_r;
public:
	void Half_Node_Flux();
};
void Half_Node_Flux();