#include <vector>
#include "2D_Euler_Solver.h"

using namespace std;

int main(int argc, char ** argv )
{
	Simulation * two_dim_Euler_Solver = new Simulation();

	two_dim_Euler_Solver->Run();

	delete two_dim_Euler_Solver;

	return 0;
}

void Simulation::Run()
{
	//建立全局参数、变量、流场
	Init_Global_Param();

	Generate_Mesh();

	Flow_Initialization();


	for (current_step = 0; current_step < max_num_of_steps; ++current_step)
	{
		//先计算x方向
		//Load_Q();

		//Half_Node_Q();

		//Half_Node_Flux();

		//Spatial_Derivative();

		//Time_Integral();
	}

	//再计算y方向

}

void Init_Global_Param() 
{
	max_num_of_steps = 10000;

	cfl_num = 0.1;

	time_step = 0.0; //时间要根据最大特征值确定，这里还没法算，只是初始化

	grid_point_num_x = 101;

	grid_point_num_y = 51;
	
	num_of_prim_vars = 4;
}

void Generate_Mesh() 
{
	double hx = 4.0, hy = 2.0;
	dx = hx / (grid_point_num_x - 1);
	dy = hy / (grid_point_num_y - 1);

	x_coord.resize(grid_point_num_x);
	y_coord.resize(grid_point_num_y);

	for (int iNode = 0; iNode < grid_point_num_x; ++iNode)
	{
		x_coord[iNode] = dx * iNode;
	}

	for (int jNode = 0; jNode < grid_point_num_y; ++jNode)
	{
		y_coord[jNode] = dy * jNode;
	}
}

void Flow_Initialization() 
{
	//流场初始化
	qField.resize(num_of_prim_vars);
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		qField[iVar].resize(grid_point_num_x);
		for (int i = 0; i < grid_point_num_x; i++)
		{
			qField[iVar][i].resize(grid_point_num_y);
		}
	}

	//标记网格点是否在流场内
	marker.resize(grid_point_num_x); 
	for (int i = 0; i < grid_point_num_x; i++)
	{
		marker[i].resize(grid_point_num_y);
	}

	for (int i = 0; i < grid_point_num_x; i++)
	{
		for (int j = 0; j < grid_point_num_y; j++)
		{
			marker[i][j] = 1;

			double x_node = x_coord[i];
			double y_node = y_coord[j];

			if (x_node > 0.6 && y_node > 0.8 && y_node < 1.2)
			{
				marker[i][j] = 0;		//标记哪些点是在物体里面，不参加计算
			}
		}
	}

	//流场赋初值
	for (int i = 0; i < grid_point_num_x; i++)
	{
		for (int j = 0; j < grid_point_num_y; j++)
		{
			if (marker[i][j] == 0) continue;

			qField[IR][i][j] = 1.0;
			qField[IU][i][j] = 3.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = 0.71429;
		}
	}
}