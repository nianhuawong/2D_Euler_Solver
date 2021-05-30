#pragma once
#include <vector>
using namespace std;

class Structured_Mesh;
extern Structured_Mesh* mesh;
extern int num_ghost_point;
extern int num_grid_point_x, num_grid_point_y;
extern int total_points_x, total_points_y;
extern int num_half_point_x, num_half_point_y;
extern int Iw, Jw1, Jw2;
extern double dx, dy;

class Point
{
public:
	Point();
	~Point() {}

protected:
	double xPoint;
	double yPoint;

public:
	void Get_Point_Coord(double* xNode, double* yNode) { *xNode = this->xPoint; *yNode = this->yPoint; }
	void Set_Point_Coord(double xNode, double  yNode) { this->xPoint = xNode;  this->yPoint = yNode; }
};

class Structured_Mesh
{
public:
	Structured_Mesh() { NI = 0; NJ = 0; }
	Structured_Mesh(int NI, int NJ);
	~Structured_Mesh() {}

protected:
	int NI, NJ;
	vector< vector< Point > > grid_points;
	vector< vector< int > > marker;

public:
	vector< vector< Point > >& Get_Grid_Points() { return grid_points; }
	vector< vector< int > >& Get_Marker() { return marker; }
	void Set_Mesh_Dimension(int NI, int NJ) { this->NI = NI; this->NJ = NJ; }
};

void Set_Mesh_Dimension(int NI, int NJ);
