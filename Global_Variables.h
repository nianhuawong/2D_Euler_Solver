#pragma once
#include<vector>
using namespace std;

namespace GLOBAL
{
	int num_of_prim_vars;

	int current_step, max_num_of_steps;

	int num_grid_point_x, num_grid_point_y, num_ghost_point;

	int total_points_x, total_points_y;

	int num_half_point_x, num_half_point_y;

	int Iw, Jw1, Jw2;

	double dx, dy;

	double cfl_num, time_step;

	int method_of_half_q;
	int method_of_limiter;
	int method_of_flux;
	double muscl_k;
	double entropy_fix_coeff;

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

template < typename T >
void Allocate_2D_Vector(vector< vector< T > > & array, int dim1, int dim2)
{
	array.resize(dim1);
	for (int i = 0; i < dim1; i++)
	{
		array[i].resize(dim2);
	}
}
}

class Point
{
public:
	Point(){}
	~Point() {}

protected:
	double xPoint;
	double yPoint;

public:
	void Get_Point_Coord(double* xNode, double* yNode) { *xNode = this->xPoint; *yNode = this->yPoint; }
	void Set_Point_Coord(double xNode,  double  yNode) { this->xPoint = xNode;  this->yPoint = yNode; }
};

class Structured_Mesh
{
public:
	Structured_Mesh();
	~Structured_Mesh() {};

protected:
	vector< vector< Point > > grid_points;
	vector< vector< int > > marker;

public:
	vector< vector< Point > >& Get_Grid_Points() { return grid_points; }
	vector< vector< int > >& Get_Marker() { return marker; }
};

Structured_Mesh* mesh;