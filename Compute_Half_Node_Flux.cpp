#include <iostream>
#include "Global_Variables.h"
#include "Compute_Half_Node_Flux.h"
using namespace GLOBAL;

void Half_Node_Flux()
{
	auto* half_node_flux = new Half_Node_Flux_Solver();

	half_node_flux->Half_Node_Flux();

	delete half_node_flux;
}

Half_Node_Flux_Solver::Half_Node_Flux_Solver()
{
	this->method_of_flux = 1;		//1-Roe, 2-WENO, 3-WCNS

	GLOBAL::half_node_flux.resize(num_of_prim_vars);
	this->half_node_flux_l.resize(num_of_prim_vars);
	this->half_node_flux_r.resize(num_of_prim_vars);

	this->M_Dim = grid_point_num_x - 1;   //半点个数，即单元数，点数减1
	this->N_Dim = grid_point_num_y - 1;
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		half_node_flux.resize		 (M_Dim);
		half_node_flux_l[iVar].resize(M_Dim);
		half_node_flux_r[iVar].resize(M_Dim);
		for (int i = 0; i < M_Dim; i++)
		{
			half_node_flux  [iVar][i].resize(N_Dim);
			half_node_flux_l[iVar][i].resize(N_Dim);
			half_node_flux_r[iVar][i].resize(N_Dim);
		}
	}
}

void Half_Node_Flux_Solver::Half_Node_Flux()
{

}