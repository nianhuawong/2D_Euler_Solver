#pragma once
#include<vector>
using namespace std;

namespace GLOBAL
{
	int num_of_prim_vars;

	int current_step, max_num_of_steps;

	int grid_point_num_x, grid_point_num_y, ghost_point_num;

	int total_points_x, total_points_y;

	int Iw, Jw1, Jw2;

	double dx, dy;

	double cfl_num, time_step;

	vector<double> x_coord, y_coord;

	vector< vector< int > > marker;
	vector< vector< vector< double > > > qField;
	vector< vector< vector< double > > > half_node_Q_l;
	vector< vector< vector< double > > > half_node_Q_r;
	vector< vector< vector< double > > > half_node_flux;

#define IR 0
#define IU 1
#define IV 2
#define IP 3

void Init_Global_Param();
void Generate_Mesh();
void Flow_Initialization();
void Compute_Boundary();

}