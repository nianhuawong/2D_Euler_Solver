#include <vector>
#include <iostream>
#include "Compute_Half_Node_Flux.h"
#include "Global_Variables.h"

void Half_Node_Flux()
{
	auto* half_node_flux = new Half_Node_Flux_Solver();

	half_node_flux->Half_Node_Flux();

	delete half_node_flux;
}

Half_Node_Flux_Solver::Half_Node_Flux_Solver()
{
	half_node_flux_l.resize(num_of_prim_vars);
	half_node_flux_r.resize(num_of_prim_vars);

	int Mdim = grid_point_num_x - 1;   //半点个数，即单元数，点数减1
	int Ndim = grid_point_num_y - 1;
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		half_node_flux_l[iVar].resize(Mdim);
		half_node_flux_r[iVar].resize(Mdim);
		for (int i = 0; i < Mdim; i++)
		{
			half_node_flux_l[iVar][i].resize(Ndim);
			half_node_flux_r[iVar][i].resize(Ndim);
		}
	}
}

void Half_Node_Flux_Solver::Half_Node_Flux()
{

}