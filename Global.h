#pragma once
#include<vector>
using namespace std;

namespace GLOBAL
{
#define IR 0
#define IU 1
#define IV 2
#define IP 3

	int num_of_prim_vars;

	int current_step, max_num_of_steps;

	double cfl_num, time_step;

	int method_of_half_q;
	int method_of_limiter;
	int method_of_flux;
	double muscl_k;
	double entropy_fix_coeff;

	template < typename T >
	void Allocate_2D_Vector(vector< vector< T > >& array, int dim1, int dim2)
	{
		array.resize(dim1);
		for (int i = 0; i < dim1; i++)
		{
			array[i].resize(dim2);
		}
	}
}
