#include <iostream>
#include "2D_Euler_Solver.h"
#include "QlQr_Solver.h"
#include "Global.h"
#include "Geometry.h"
#include <fstream>

#ifndef _WIN32
#include <sys/time.h>
#include <unistd.h>
#include <sys/times.h>
#endif

#ifdef WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

int			num_of_prim_vars;
int			current_step, max_num_of_steps;
double		cfl_num, time_step, physical_time, max_simu_time;
int			method_of_half_q;
int			method_of_limiter;
int			method_of_flux;
double		muscl_k;
double		entropy_fix_coeff;
char		solve_direction;
int			residual_output_steps;
int			flow_save_steps;
double		converge_criterion;
string		tec_file_name;
int			num_of_RK_stages;
VDouble2D	RK_Coeff;
clock_t		lastTime, nowTime;

void Init_Global_Param()
{
	//==============================================================================================
	lastTime = Get_Current_Time(); 
	nowTime  = lastTime;

	num_of_prim_vars	= 4;			//原始变量个数，控制方程个数
	physical_time		= 0.0;
	time_step			= 0.0;			//时间步长要根据最大特征值确定，这里只是初始化
	solve_direction		= 'x';
	
	num_of_RK_stages	= 3;
	RK_Coeff			= { {1.0, 0.0, 1.0},{3.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0},{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0} };
	
	//num_of_RK_stages	= 2;
	//RK_Coeff			= { {1.0, 0.0, 1.0},{0.5, 0.5, 0.5} };
	//==============================================================================================

	num_grid_point_x = 2 * 240 + 1;
	num_grid_point_y = 2 * 60  + 1;

	max_num_of_steps			= 2000;
	residual_output_steps		= 2;		//残差输出间隔步数
	flow_save_steps				= 50;		//流场输出间隔步数
	converge_criterion			= 1e-8;		//残差收敛标准
	tec_file_name				= "./results/flow.plt";

	cfl_num						= 0.3;
	max_simu_time				= 0.2;

	method_of_half_q			= 3;		//1-MUSCL,	  2-WENO(不插值),   3-WCNS
	muscl_k						= 1.0/3;	//0.0-二阶迎风偏置，		    1/3-二阶迎风偏置
	method_of_limiter			= 1;		//0-nolim,    1-vanleer,        2-minmod,	  3-superbee,	4-1st
	method_of_flux				= 2;		//1-Roe,	  2-Steger,			3-VanLeer,    4-WENO,		5-WCNS 
	entropy_fix_coeff			= 0.001;	//Roe格式熵修正系数epsilon
	//==============================================================================================
	MakeDirectory("./results");

#ifndef _WIN32
	Read_Parameter_File("./input.txt");	//在集群上计算时，可通过参数文件设置参数
#endif // !_WIN32	
}

void Load_Q()
{
	qField = qField_N1;
}

void Set_Solve_Direction(char direction)
{
	solve_direction = direction;
	if (direction=='x')
	{		
		qField_N0 = qField;		//RK公式里第一项，当前时间步Q
		qField_N1 = qField;		//RK公式里第二项，下一stage的Q
	}
	else if(direction == 'y')
	{
		//最关键的点：y方向计算的qField和x方向计算的qField是相同的！！这是算子分裂法的关键。
		//也正因为如此，在计算y方向时，无需重新计算边界条件
		qField	  = qField_N0;		//将qField还原为计算x方向之前的Q值，y方向还是用原来的qField来计算rhs
		qField_N0 = qField_N1;		//RK公式里第一项，是x方向多个stage推进完后，求出来的Q值
		qField_N1 = qField;			//qField和qField_N1是RK推进中的关键变量，仍然设为计算x方向之前的Q值
	}
}

void Primitive_To_Conservative( VDouble & primitive, VDouble &conservative )
{
	double gama = 1.4;
	double rho, u, v, p;
	rho = primitive[IR];
	u = primitive[IU];
	v = primitive[IV];
	p = primitive[IP];

	conservative[IR] = rho;
	conservative[IU] = rho * u;
	conservative[IV] = rho * v;
	conservative[IP] = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
}

void Conservative_To_Primitive(VDouble& conservative, VDouble &primitive)
{
	double gama = 1.4;
	double rho		= conservative[IR];
	double rho_u	= conservative[IU];
	double rho_v	= conservative[IV];
	double E		= conservative[IP];

	double u = rho_u / rho;
	double v = rho_v / rho;
	double p = (gama - 1) * (E - 0.5 * rho * (u * u + v * v));

	primitive[IR] = rho;
	primitive[IU] = u;
	primitive[IV] = v;
	primitive[IP] = p;
}

//计算内场点的标号范围
void Get_IJK_Region(int& ist, int& ied, int& jst, int& jed)
{
	ist = num_ghost_point;
	ied = num_ghost_point + num_grid_point_x - 1;//ist->ied,范围包含了ied

	jst = num_ghost_point;
	jed = num_ghost_point + num_grid_point_y - 1;//jst->jed,范围包含了jed
}

bool Need_Stop_Iteration()
{
	return stop_by_residual || physical_time >= max_simu_time || current_step == max_num_of_steps;
}

bool IsNaN(VDouble& data)
{
	bool flag = 0;
	for (int i = 0; i < data.size(); i++)
	{
		if (data[i] != data[i])
		{
			flag = 1;
			break;
		}
	}
	return flag;
}

clock_t Get_Current_Time()
{
#ifdef _WIN32
	return clock();
#else
	struct tms tp;
	times(&tp);
	return tp.tms_stime;	//system time
	//return tp.tms_utime;	//cpu time
#endif
}

double GetClockTicksPerSecond()
{
	long clockTicksPerSecond;

#ifdef _WIN32
	clockTicksPerSecond = CLOCKS_PER_SEC;
#else
	clockTicksPerSecond = sysconf(_SC_CLK_TCK);
#endif
	return clockTicksPerSecond;
}

void ExtractValue(VDouble primitiveVector, double& rm, double& um, double& vm, double& pm)
{
	rm = primitiveVector[IR];
	um = primitiveVector[IU];
	vm = primitiveVector[IV];
	pm = primitiveVector[IP];
}

void Read_Parameter_File(string fileName)
{
	fstream file;
	file.open(fileName, ios_base::in);

	int n;
	file >> n;
	num_grid_point_x = n * 240 + 1;
	num_grid_point_y = n * 60  + 1;

	file >> method_of_half_q;	//半节点插值方法	//1-MUSCL,	  2-WENO(不插值),   3-WCNS
	file >> method_of_limiter;	//MUSCL限制器种类	//0-nolim,    1-vanleer,        2-minmod,	  3-superbee,	4-1st
	file >> method_of_flux;		//通量计算方法		//1-Roe,	  2-Steger,			3-VanLeer,    4-WENO,		5-WCNS

	file.clear();
	file.close();
}

void MakeDirectory(const string& directoryName)
{
#ifdef WIN32
	int flag = _mkdir(directoryName.c_str());
#else
	int flag = mkdir(directoryName.c_str(), S_IRWXU);
#endif
	if (flag == 0)
	{
		cout << directoryName << " directory has been created successfully !\n";
	}
}