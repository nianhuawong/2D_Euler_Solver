#pragma once

#include <iostream>

class Simulation
{
public:
	Simulation() {}
	~Simulation(){}

public:
	void Run();

};

int current_step,		max_num_of_steps;
int grid_point_num_x,	grid_point_num_y, ghost_point_num;
int total_points_x,		total_points_y;
int Iw, Jw1, Jw2;
int num_of_prim_vars;
double dx, dy;
double cfl_num, time_step;
vector<double> x_coord, y_coord;
vector< vector< vector< double > > > qField;
vector< vector< int > > marker;

#define IR 0
#define IU 1
#define IV 2
#define IP 3

void Init_Global_Param();
void Generate_Mesh();
void Flow_Initialization();
void Compute_Boundary();