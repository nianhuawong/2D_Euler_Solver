#include "Post_Process.h"
#include <fstream>
#include <iostream>
#include "Global.h"
#include "QlQr_Solver.h"
#include "Geometry.h"
#include <iomanip>

bool stop_by_residual = 0;
void Compute_Residual()
{
	auto* residual_solver = new Residual();

	residual_solver->Compute_Residual();

	delete residual_solver;
}

Residual::Residual()
{
	res_L1.resize (num_of_prim_vars);
	res_L2.resize (num_of_prim_vars);
	res_Loo.resize(num_of_prim_vars);
}

void Residual::Compute_Residual()
{
	vector< vector< int > >& marker = mesh->Get_Marker();
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		int count = 0;
		double res_max = -1e40;
		vector< vector< double > >& qv = qField[iVar];
		for (int i = 0; i < num_half_point_x; i++)
		{
			for (int j = 0; j < num_half_point_y; j++)
			{
				if (marker[i][j] == 0) continue;
				
				double res1 = abs( qField_N1[iVar][i][j] - qField[iVar][i][j] );

				double res2 = res1 * res1;

				res_max = res1 > res_max ? res1 : res_max;

				res_L1[iVar] += res1;
				res_L2[iVar] += res2;
				res_Loo[iVar] = res_max;

				count++;
			}
		}
		res_L1[iVar] /= count;
		res_L2[iVar] = sqrt(res_L2[iVar] / count);
	}
	
	this->OutputResidual();
}
void Residual::OutputResidual()
{
	bool flag = (current_step + 1) % residual_output_steps;//整除时flag=0,输出残差，非整除时flag=1，不输出
	if (flag) return;

	cout << setiosflags(ios::right);
	cout << setiosflags(ios::scientific);
	cout << setprecision(5);
	cout << "Iteration        rho_res_Loo      u_res_Loo     v_res_Loo    p_res_Loo" << endl;
	cout << current_step + 1 << "    "
		 << res_Loo[IR]		 << "    " << res_Loo[IU] << "    "
		 << res_Loo[IV]		 << "    " << res_Loo[IP] << endl;

	if (res_Loo[IR] < converge_criterion &&
		res_Loo[IU] < converge_criterion &&
		res_Loo[IV] < converge_criterion &&
		res_Loo[IP] < converge_criterion   )
	{
		stop_by_residual = 1;
	}
}

bool Residual::Stop_by_Residual()
{
	return false;
}

bool Stop_by_Residual()
{
	return false;
}

void Output_Flowfield()
{
	bool flag = (current_step + 1) % flow_save_steps;//整除时flag=0,输出，非整除时flag=1，不输出
	if (flag) return;

	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();
	vector< vector< int > >& marker = mesh->Get_Marker();

	cout << "dumping flowfield..." << endl;
	fstream file;
	file.open("./tecflow/flow.plt", ios_base::out);

	file << "VARIABLES         = \"x\", \"y\", \"rho\", \"u\",\"v\", \"p\"" << endl;
	file << "ZONE T = \"Zone 1\"" << endl;
	file << "I = " << num_grid_point_x << "  J = " << num_grid_point_y << "  K =1" << "  ZONETYPE=Ordered" << endl;
	file << "DATAPACKING=POINT"   << endl;

	file << setiosflags(ios::right);
	file << setiosflags(ios::scientific);
	file << setprecision(8);

	for (int i = 0; i < num_grid_point_x; i++)
	{
		for (int j = 0; j < num_grid_point_y; j++)
		{
			double xcoord, ycoord;
			grid_points[i][j].Get_Point_Coord(xcoord, ycoord);

			double rho = qField_N1[IR][i][j];
			double u   = qField_N1[IU][i][j];
			double v   = qField_N1[IV][i][j];
			double E   = qField_N1[IP][i][j];
			double p   = Energy_2_Pressure(E, rho, u, v);

			file << xcoord << "    " << ycoord << "    " << rho << "    " << u << "    " 
				 << v      << "    " << p << endl;
		}
	}
	file.close();
}