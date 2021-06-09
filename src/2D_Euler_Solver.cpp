//////////////////////////////////////////////////////////
//														//
//			2D Euler Solver on Structured Grid			//
//					by Wang Nianhua						//
//				Email: nianhuawong@126.com				//
//					   2021.6.7							//
//////////////////////////////////////////////////////////
#include "2D_Euler_Solver.h"
#include "Global.h"

#ifdef _OPENMP
#include <iostream>
#include <omp.h>
#include <cstdlib>
#endif

int main(int argc, char ** argv )
{
#ifdef _OPENMP
	InitializeOpenMP(argv);
#endif
	Simulation * two_dim_Euler_Solver = new Simulation();

	two_dim_Euler_Solver->Run();

	delete two_dim_Euler_Solver;

	return 0;
}

void Simulation::Run()
{
	Init_Global_Param();

	Generate_Mesh();

	Init_Flow();

	for (current_step = 1; current_step <= max_num_of_steps; ++current_step)
	{
		Solve_Time_Step();

		Compute_Boundary();

		Set_Solve_Direction('x');
		Time_Integration();

		Set_Solve_Direction('y');
		Time_Integration();

		Post_Solve();

		if (Need_Stop_Iteration())
		{
			break;
		}

		Set_Field();
	}

	Output_Flowfield();	
}

#ifdef _OPENMP
void InitializeOpenMP(char** argv)
{
	int num_threads = 1;//无命令行参数，默认为单线程
	if (argv[1] != NULL)
	{
		num_threads = atoi(argv[1]);
	}

#ifdef WIN32	
	omp_set_num_threads(num_threads);
#else	
	int omp_num_threads = omp_get_num_threads();//环境变量OMP_NUM_THREADS
	if (omp_num_threads != num_threads && argv[1] != NULL)
	{
		cout << "同时存在环境变量和命令行参数，以命令行线程数为准！" << endl;
		omp_set_num_threads(num_threads);
	}
#endif
#pragma omp parallel 
	{
		if (omp_get_thread_num() == 0)
		{
			cout << "OpenMP线程级并行已开启，线程总数为: num_of_threads =" 
				 << omp_get_num_threads() << endl;
		}
		//cout << "Hello world from open mp thread " << omp_get_thread_num() << endl;
	}
}
#endif
