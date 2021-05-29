#include "Geometry.h"
#include "Global.h"
#include "2D_Euler_Solver.h"
using namespace GLOBAL;

Point::Point()
{
	xPoint = 0.0;
	yPoint = 0.0;
}

void GLOBAL::Generate_Mesh()
{
	mesh = new Structured_Mesh(num_grid_point_x, num_grid_point_y);
	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();
	vector< vector< int > >& marker = mesh->Get_Marker();

	double hx = 4.0, hy = 2.0;
	double hx1 = 0.6, hy1 = 0.8, hy2 = 1.2;

	num_ghost_point = 2;

	total_points_x = num_grid_point_x + 2 * num_ghost_point;
	total_points_y = num_grid_point_y + 2 * num_ghost_point;

	num_half_point_x = num_grid_point_x - 1;   //半点个数，即单元数，点数减1
	num_half_point_y = num_grid_point_y - 1;   //半点个数，即单元数，点数减1
	
	Iw  = hx1 / hx * (num_grid_point_x - 1); //物面左边界的标号I，起始标号为0，不包含虚拟点
	Jw1 = hy1 / hy * (num_grid_point_y - 1); //物面下边界的标号J
	Jw2 = hy2 / hy * (num_grid_point_y - 1); //物面上边界的标号J

	dx = hx / (num_grid_point_x - 1);
	dy = hy / (num_grid_point_y - 1);

	for (int i = 0; i < num_grid_point_x; i++)
	{
		for (int j = 0; j < num_grid_point_y; j++)
		{
			double x_node = i * dx;
			double y_node = j * dy;

			grid_points[i][j].Set_Point_Coord(x_node, y_node);

			marker[i][j] = 1;
			if (x_node > hx1 && y_node > hy1 && y_node < hy2)
			{
				marker[i][j] = 0;		//标记哪些点是在物体里面，不参加计算
			}
		}
	}
}

Structured_Mesh::Structured_Mesh(int NI, int NJ)
{
	this->NI = NI;
	this->NJ = NJ;
	Allocate_2D_Vector(grid_points, NI, NJ);
	Allocate_2D_Vector(marker, NI, NJ);
}
