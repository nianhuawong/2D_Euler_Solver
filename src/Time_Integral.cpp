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
	//准备开始写3阶TVD RK方法
	//
	Solve_QlQr();

	Solve_Flux();

	Solve_Spatial_Derivative();
	//

	Solve_Time_Step();

	Update_Flowfield();
}

void Update_Flowfield()
{
	//if (solve_direction == 'x')
	//{
		//Update_Flowfield_X();
	//}
	//else if (solve_direction == 'y')
	//{
		Update_Flowfield_Y();
	//}
}
void Update_Flowfield_X()
{
	VInt2D& marker = mesh->Get_Marker();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	VDouble qPrimitive(num_of_prim_vars);
	VDouble qConservative(num_of_prim_vars);
	VDouble rhsVector(num_of_prim_vars);
	for (int j = jst; j < jed; j++)
	{
		for (int i = ist; i < ied; i++)
		{
			if (marker[i][j] == 0) continue;

			qPrimitive = qField[i][j];
			rhsVector = rhs[i][j];
			Primitive_To_Conservative(qPrimitive, qConservative);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				qConservative[iVar] += time_step * rhsVector[iVar];
			}

			Conservative_To_Primitive(qConservative, qPrimitive);
			qField_N1[i][j] = qPrimitive;
		}
	}
}

void Update_Flowfield_Y()
{
	VInt2D& marker = mesh->Get_Marker();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	VDouble qPrimitive   (num_of_prim_vars);
	VDouble qConservative(num_of_prim_vars);
	VDouble rhsVector    (num_of_prim_vars);
	for (int j = jst; j < jed; j++)
	{
		for (int i = ist; i < ied; i++)
		{
			if (marker[i][j] == 0) continue;

			qPrimitive = qField_N1[i][j];
			rhsVector  = rhs[i][j];
			Primitive_To_Conservative(qPrimitive, qConservative);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				qConservative[iVar] += time_step * rhsVector[iVar];
			}

			Conservative_To_Primitive(qConservative, qPrimitive);
			qField_N1[i][j] = qPrimitive;
			if (qPrimitive[IR] != qPrimitive[IR])
			{
				int kkk = 1;
			}
		}
	}

	//qField_N2[iVar][i][j] = 3.0 / 4.0 * qField[iVar][i][j] + 1.0 / 4.0 * qField_N1[iVar][i][j] + 1.0 / 4.0 * time_step * rhs1[iNode];

	//qField_N3[iVar][i][j] = 1.0 / 3.0 * qField[iVar][i][j] + 2.0 / 3.0 * qField_N2[iVar][i][j] + 2.0 / 3.0 * time_step * rhs2[iNode];

}