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

	Get_IJK_Region(ist, ied, jst, jed);
}

void Residual::Compute_Residual()
{
	VInt2D& marker = mesh->Get_Marker();
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		int count = 0;
		double res_max = -1e40;
		
		for (int i = ist; i < ied; i++)
		{
			for (int j = jst; j < jed; j++)
			{
				if (marker[i][j] == 0) continue;
				
				double res1 = abs( qField_N1[i][j][iVar] - qField[i][j] [iVar]);

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
	bool flag1 = current_step % residual_output_steps;//整除时flag=0,输出残差，非整除时flag=1，不输出
	if (flag1) return;

	cout << setiosflags(ios::left);
	cout << setiosflags(ios::scientific);
	cout << setprecision(5);

	bool flag2 = current_step % (10*residual_output_steps);
	if (!flag2)
	{
		cout << "Iteration\trho_res_Loo\tu_res_Loo\tv_res_Loo\tp_res_Loo" << endl;
	}
	
	cout << current_step + 1 << "\t      "
		 << res_Loo[IR]		 << "\t" << res_Loo[IU] << "\t"
		 << res_Loo[IV]		 << "\t" << res_Loo[IP] << endl;

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
	bool flag = current_step % flow_save_steps;//整除时flag=0,输出，非整除时flag=1，不输出
	if (flag) return;

	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();
	VInt2D& marker = mesh->Get_Marker();

	cout << "dumping flowfield..." << endl;
	fstream file;
	file.open("./flow.plt", ios_base::out);

	file << "VARIABLES = \"x\", \"y\", \"rho\", \"u\",\"v\", \"p\"" << endl;
	file << "ZONE T = \"Zone 1\"" << endl;
	file << "I = " << num_grid_point_x << "  J = " << num_grid_point_y << "  K =1" << "  ZONETYPE=Ordered" << endl;
	file << "DATAPACKING=POINT"   << endl;

	file << setiosflags(ios::right);
	file << setiosflags(ios::scientific);
	file << setprecision(8);

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	for (int j = jst; j < jed; j++)
	{
		for (int i = ist; i < ied; i++)
		{
			//if (marker[i][j] == 0) continue;

			double xcoord, ycoord;
			grid_points[i][j].Get_Point_Coord(xcoord, ycoord);

			double rho = qField_N1[i][j][IR];
			double u   = qField_N1[i][j][IU];
			double v   = qField_N1[i][j][IV];
			double p   = qField_N1[i][j][IP];

			file << xcoord << "    " << ycoord << "    " << rho << "    " << u << "    " 
				 << v      << "    " << p << endl;
		}
	}
	
	file.close();
}