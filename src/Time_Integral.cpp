#include "2D_Euler_Solver.h"
#include "Spatial_Derivative.h"
#include "Geometry.h"
#include "Global.h"
#include "Flux_Solver.h"
#include "QlQr_Solver.h"
#include "Time_Integral.h"


void Time_Integration()
{
	auto* time_marching = new Time_Marching_Solver();
	time_marching->Time_Marching();
	delete time_marching;
}

Time_Marching_Solver::Time_Marching_Solver()
{

}

void Time_Marching_Solver::Time_Marching()
{
	Solve_QlQr();

	Solve_Flux();

	Solve_Spatial_Derivative();

	Solve_Time_Step();

	Update_Flowfield();
}

void Update_Flowfield()
{
	if (solve_direction == 'x')
	{
		Update_Flowfield_X();
	}
	else if (solve_direction == 'y')
	{
		Update_Flowfield_Y();
	}
}
void Update_Flowfield_X()
{
	VInt2D& marker = mesh->Get_Marker();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	VDouble qPrimitive(num_of_prim_vars);
	VDouble qConservative(num_of_prim_vars);
	
	for (int j = jst; j < jed; j++)
	{
		for (int i = ist; i < ied; i++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				qPrimitive = { qField[IR][i][j],qField[IU][i][j],qField[IV][i][j],qField[IP][i][j] };
				Primitive_To_Conservative(qPrimitive, qConservative);

				qConservative[iVar]  += time_step * rhs[iVar][i][j];
			}

			Conservative_To_Primitive(qConservative, qPrimitive);
			qField_N1[IR][i][j] = qPrimitive[IR];
			qField_N1[IU][i][j] = qPrimitive[IU];
			qField_N1[IV][i][j] = qPrimitive[IV];
			qField_N1[IP][i][j] = qPrimitive[IP];
		}
	}

//qField_N2[iVar][i][j] = 3.0 / 4.0 * qField[iVar][i][j] + 1.0 / 4.0 * qField_N1[iVar][i][j] + 1.0 / 4.0 * time_step * rhs1[iNode];

//qField_N3[iVar][i][j] = 1.0 / 3.0 * qField[iVar][i][j] + 2.0 / 3.0 * qField_N2[iVar][i][j] + 2.0 / 3.0 * time_step * rhs2[iNode];

}

void Update_Flowfield_Y()
{
	VInt2D& marker = mesh->Get_Marker();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	VDouble qPrimitive(num_of_prim_vars);
	VDouble qConservative(num_of_prim_vars);

	for (int j = jst; j < jed; j++)
	{
		for (int i = ist; i < ied; i++)
		{
			if (marker[i][j] == 0) continue;
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				qPrimitive = { qField_N1[IR][i][j],qField_N1[IU][i][j],qField_N1[IV][i][j],qField_N1[IP][i][j] };
				Primitive_To_Conservative(qPrimitive, qConservative);

				qConservative[iVar] += time_step * rhs[iVar][i][j];
			}

			Conservative_To_Primitive(qConservative, qPrimitive);
			qField_N1[IR][i][j] = qPrimitive[IR];
			qField_N1[IU][i][j] = qPrimitive[IU];
			qField_N1[IV][i][j] = qPrimitive[IV];
			qField_N1[IP][i][j] = qPrimitive[IP];
		}
	}

	//qField_N2[iVar][i][j] = 3.0 / 4.0 * qField[iVar][i][j] + 1.0 / 4.0 * qField_N1[iVar][i][j] + 1.0 / 4.0 * time_step * rhs1[iNode];

	//qField_N3[iVar][i][j] = 1.0 / 3.0 * qField[iVar][i][j] + 2.0 / 3.0 * qField_N2[iVar][i][j] + 2.0 / 3.0 * time_step * rhs2[iNode];

}