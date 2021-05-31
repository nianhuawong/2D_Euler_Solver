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
	void Get_Point_Coord(double& xNode, double& yNode) { xNode = this->xPoint; yNode = this->yPoint; }
	void Set_Point_Coord(double xNode, double  yNode) { this->xPoint = xNode;  this->yPoint = yNode; }
};

class Structured_Mesh
{
public:
	Structured_Mesh(int NI_tot, int NJ_tot);
	~Structured_Mesh() {}

protected:
	int NI, NJ;
	int NI_tot, NJ_tot;
	int ist, ied, jst, jed;
	int num_ghost_point;
	double dx, dy;
	vector< vector< Point > > grid_points;
	VInt2D marker;

public:
	vector< vector< Point > >& Get_Grid_Points() { return grid_points; }
	VInt2D& Get_Marker() { return marker; }
	void Set_Mesh_Dimension(int NI, int NJ) { this->NI = NI; this->NJ = NJ; }
	void Set_Dxdy(double dx, double dy) { this->dx = dx; this->dy = dy; }
	void Set_Range(int ist, int ied, int jst, int jed);
	void Set_Num_Ghost_Point(int num_ghost_point) { this->num_ghost_point = num_ghost_point; }
};

void Set_Mesh_Dimension(int NI, int NJ);
void Get_IJK_Region(int& ist, int& ied, int& jst, int& jed);