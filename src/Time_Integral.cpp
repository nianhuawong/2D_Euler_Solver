#include "2D_Euler_Solver.h"
#include "Spatial_Derivative.h"
#include "Geometry.h"
#include "Global.h"
#include "Flux_Solver.h"
#include "QlQr_Solver.h"
#include "Time_Integral.h"


void Time_Integration()
{
	Time_Marching_Solver* time_marching = new Time_Marching_Solver();
	time_marching->Time_Marching();
	delete time_marching;
}

Time_Marching_Solver::Time_Marching_Solver()
{

}

void Time_Marching_Solver::Time_Marching()
{
	for (int iStage = 0; iStage < num_of_RK_stages; iStage++)
	{
		Load_Q();

		Solve_QlQr();

		Solve_Flux();

		Solve_Spatial_Derivative();

		Update_Flowfield(iStage);
	}	
}

void Update_Flowfield(int iStage)
{
	double RK_Coeff_a = RK_Coeff[iStage][0];
	double RK_Coeff_b = RK_Coeff[iStage][1];
	double RK_Coeff_c = RK_Coeff[iStage][2];

	VInt2D& marker = mesh->Get_Marker();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	VDouble rhsVector	  (num_of_prim_vars);
	VDouble qPrimitive0   (num_of_prim_vars);
	VDouble qPrimitive1   (num_of_prim_vars);
	VDouble qConservative0(num_of_prim_vars);
	VDouble qConservative1(num_of_prim_vars);

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif	
	for (int j = jst; j < jed; j++)
	{
		for (int i = ist; i < ied; i++)
		{
			if (marker[i][j] == 0) continue;

			rhsVector = rhs[i][j];			//rhs(q0),用上一步的q计算得到的rhs

			qPrimitive0 = qField_N0[i][j];	//RK公式里第一项,q0
			Primitive_To_Conservative(qPrimitive0, qConservative0);

			qPrimitive1 = qField_N1[i][j];	//RK公式里第二项，q0、q1、q2，也即上一stage的q值
			Primitive_To_Conservative(qPrimitive1, qConservative1);
			
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				qConservative1[iVar] = RK_Coeff_a * qConservative0[iVar] + RK_Coeff_b * qConservative1[iVar] + RK_Coeff_c * time_step * rhsVector[iVar];
			} 

			Conservative_To_Primitive(qConservative1, qPrimitive1);
			
			//qPrimitive1[IP] = fabs(qPrimitive1[IP]);	//保证压力不出负
			//qPrimitive1[IP] = min(qPrimitive1[IP], 5.0* qPrimitive0[IP]);

			//RK公式里左端项，q1、q2、q3，即下一stage的q值，还要继续用该值计算rhs(q1)、rhs(q2)
			qField_N1[i][j] = qPrimitive1;	
			if ( IsNaN(qPrimitive1) )
			{
				int kkk = 1;
			}
		}
	}
}

void Set_Field()
{
	//多步RK推进完成之后，流场变量更新在qField_N1中，要重新返回qField，以便下一步时间步迭代
	qField = qField_N1;
}