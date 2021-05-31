#pragma once
#include<string>
class Simulation
{
public:
	Simulation() {}
	~Simulation(){}

public:
	void Run();
};

void Init_Global_Param();
void Generate_Mesh();
void Init_Flow();
void Compute_Boundary();
void Load_Q();
void Set_Solve_Direction(char direction);
void Solve_QlQr();
void Solve_Flux();
void Solve_Spatial_Derivative();
void Solve_Time_Step();
void Time_Integration();
void Compute_Residual();
bool Stop_by_Residual();
void Output_Flowfield();

void Test();