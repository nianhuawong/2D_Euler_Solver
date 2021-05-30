#pragma once

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
void Flow_Init();
void Compute_Boundary();
void Load_Q();
void Solve_QlQr();
void Solve_Flux();
void Spatial_Derivative();
void Time_Integral();

void Test();