#include <fstream>
#include <iostream>
#include <cmath>
#include "Post_Process.h"
#include "Global.h"
#include "QlQr_Solver.h"
#include "Geometry.h"
#include <iomanip>
#include "Spatial_Derivative.h"

bool stop_by_residual = 0;
void Compute_Residual()
{
	Residual* residual_solver = new Residual();

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
				
				double res1 = fabs( qField_N1[i][j][iVar] - qField[i][j] [iVar]);
				//double res1 = fabs(rhs[i][j][iVar]);
				double res2 = res1 * res1;
				res_max = max(res1, res_max);

				res_L1 [iVar] += res1;
				res_L2 [iVar] += res2;
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
	bool flag1 = current_step % residual_output_steps;//flag1=0,整除时输出残差；flag1=1，非整除时不输出
	bool flag2 = current_step == 1;					  //flag2=0,中间步不输出，  flag2=1，首步输出
	if (flag1 && flag2 == 0) return;

	if (flag2)
	{
		cout << "Iteration\trho_res_L2\tu_res_L2\tv_res_L2\tp_res_L2" << endl;
	}
	cout << setiosflags(ios::left);
	cout << setiosflags(ios::scientific);
	cout << setprecision(5);
	cout << current_step	<< "\t        "
		 << res_L2[IR]		<< "\t" << res_L2[IU] << "\t"
		 << res_L2[IV]		<< "\t" << res_L2[IP] << endl;

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
	bool flag0 = current_step % flow_save_steps;//flag0=0,整除时，输出；		  flag0=1，非整除时，不输出
	bool flag1 = Need_Stop_Iteration();			//flag1=1,退出迭代，需要输出流场；flag1=0,中间步，不输出
	if (flag0 && flag1==0) return;

	cout << "dumping flowfield..."  << "\tIter = " << current_step << "\tphysical_time =" << physical_time << endl << endl;
	if (flag1 == 0) //flag1=0,中间步，输出抬头
	{
		cout << "Iteration\trho_res_Loo\tu_res_Loo\tv_res_Loo\tp_res_Loo" << endl;
	}
	
	fstream file;
	file.open(tec_file_name, ios_base::out);

	file << "VARIABLES = \"x\", \"y\", \"rho\", \"u\",\"v\", \"p\"" << endl;
	file << "ZONE T = \"Zone 1\"" << endl;
	file << "I = " << num_grid_point_x << "  J = " << num_grid_point_y << "  K =1" << "  ZONETYPE=Ordered" << endl;
	file << "DATAPACKING=POINT"   << endl;

	file << setiosflags(ios::right);
	file << setiosflags(ios::scientific);
	file << setprecision(8);

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();
	VInt2D& marker = mesh->Get_Marker();
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

void Post_Solve()
{
	Compute_Residual();

	Output_Flowfield();
}