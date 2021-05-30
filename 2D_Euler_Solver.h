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
void Solve_QlQr();
void Solve_Flux();
